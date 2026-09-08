#!/usr/bin/env bash
# Run one CLAMPbase model through the GLS pipeline + store build, with SERIAL
# step 6 (n_workers_step6=1 -> in-process, no fork -> no ProcessPoolExecutor
# deadlock). Idempotent + resume-safe.
#
# Flow: `shortcut gls --unlock` registers the model and writes pipeline_config
# WITHOUT running compute; we then patch the config to serial step 6 and run
# `workflow gls run` (which honours the edited config, unlike shortcut which
# regenerates it). Finally `store build`.
#
# Usage: _run_clampbase_model.sh <name> <rds_path>
set -uo pipefail
name="$1"; rds="$2"

WS="${PHENOPLIER_HOME:-/pividori_lab/phenoplier_workspace}"
COHORT=phenomexcan_rapid_gwas
CENTRAL="${CLAMPBASE_CENTRAL:?set CLAMPBASE_CENTRAL to the study collection dir}"
NJOBS="${CLAMPBASE_NJOBS:-${SLURM_CPUS_PER_TASK:-12}}"
# Worker/BLAS pinning (env-overridable). Defaults suit a SLURM cgroup (pico/alpine):
# serial step6, step7 = 6 workers x 2 BLAS = 12 threads. No-cgroup hosts (local/lab)
# should export N_WORKERS_STEP7=12 BLAS_STEP7=1 and OMP_NUM_THREADS=1 (single-threaded
# workers) -- config blas_threads may be a no-op once numpy is imported.
NW_S6="${CLAMPBASE_NW_STEP6:-1}"; BL_S6="${CLAMPBASE_BLAS_STEP6:-8}"
NW_S7="${CLAMPBASE_NW_STEP7:-6}"; BL_S7="${CLAMPBASE_BLAS_STEP7:-2}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export PHENOPLIER_HOME="$WS" PHENOPLIER_ROOT_DIR="$WS"

proj="$WS/projects/$name"
cfg="$proj/pipeline_config.yaml"
summary="$proj/results/gls/phenoplier/gls-summary-${COHORT}.tsv.gz"
store="$CENTRAL/stores/${name}.h5"
mkdir -p "$CENTRAL/stores" "$CENTRAL/summaries"

echo "[start] $(date -Is)  $name  njobs=$NJOBS  mem=${SLURM_MEM_PER_NODE:-?}MB  ver=$(phenoplier -v 2>&1 | tail -1)"

if [ ! -s "$summary" ]; then
  # Phase 1: register model + generate config, NO compute (--unlock returns
  # after init). ONLY on first run (no config yet) -- on requeue we must NOT
  # re-init, so the expensive step1-6 sentinels are preserved and we resume.
  if [ ! -s "$cfg" ]; then
    phenoplier shortcut gls --input "$rds" --name "$name" --model-namespace "${CLAMPBASE_NS:-clampbase}" \
      --trait-filter biomedical --lv-percentile 0.01 --executor local --n-jobs "$NJOBS" --unlock
  fi

  # Phase 2: pin worker/BLAS counts so each job uses exactly its cgroup cores.
  #   step6 serial (deadlock-proof); step7 = 6 workers x 2 BLAS = 12 threads.
  # Without the step7 cap, cpu_budget over-detects on this node's cgroup and each
  # job spawns ~29 threads -> node load ~380 -> step7 thrashes (~10x slowdown).
  python - "$cfg" "$NW_S6" "$BL_S6" "$NW_S7" "$BL_S7" <<'PY'
import sys, yaml
p, nw6, bl6, nw7, bl7 = sys.argv[1], *map(int, sys.argv[2:6])
c = yaml.safe_load(open(p))
c["n_workers_step6"] = nw6; c["blas_threads_step6"] = bl6
c["n_workers_step7"] = nw7; c["blas_threads_step7"] = bl7
yaml.safe_dump(c, open(p, "w"), sort_keys=False)
print(f"patched step6={nw6}x{bl6} step7={nw7}x{bl7} ->", p)
PY

  # Phase 3: run the pipeline (resume-safe; clear any stale lock first).
  rm -f "$proj/.snakemake/locks/"*.lock 2>/dev/null || true
  phenoplier workflow gls run --config-file "$cfg"
else
  echo "[gls] skip (summary exists)"
fi

# Gate 1: the combined summary must EXIST and be COMPLETE (2366 unique phenotypes,
# no zero-p, degenerate-only NaN). summarize concatenates whatever per-phenotype
# files exist, so an interrupted step7 yields a smaller-but-valid-looking summary.
if [ ! -s "$summary" ]; then
  echo "[ERROR] $name: GLS produced no combined summary (run failed?). Refusing to build a store." >&2
  exit 1
fi
if ! python "$HERE/verify_summary.py" "$summary"; then
  echo "[ERROR] $name: summary incomplete/invalid. Refusing to build a store; requeue to finish step7." >&2
  exit 1
fi

# Gate 2: build the store atomically -- to a .tmp, verify it opens, then rename.
# A failed store build must NOT leave a partial .h5 that the -s guard would reuse.
if [ ! -s "$store" ]; then
  echo "[store] $(date -Is)"
  tmp="${store}.tmp.$$"
  rm -f "$tmp"
  if ! phenoplier store build --output "$tmp" --cohort "$COHORT" --clamp-rds "$rds" \
        --gls-dir "$proj/results/gls/phenoplier/${COHORT}" --method-name CLAMP; then
    echo "[ERROR] $name: store build failed." >&2; rm -f "$tmp"; exit 1
  fi
  if ! python -c "import h5py,sys; h5py.File(sys.argv[1],'r').close()" "$tmp"; then
    echo "[ERROR] $name: built store failed to open." >&2; rm -f "$tmp"; exit 1
  fi
  mv -f "$tmp" "$store"
else
  echo "[store] skip (exists)"
fi

# Normalize the (verified) summary into the central collection alongside the store.
cp -f "$summary" "$CENTRAL/summaries/${name}.tsv.gz"
echo "[done] $(date -Is)  $name"
