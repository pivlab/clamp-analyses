#!/usr/bin/env bash
# Run pivlab/phenoplier-cli's GLS pipeline for ONE CLAMP model and publish its
# combined summary where this repo's rules expect it. Driven by
# workflow/rules/phenoplier.smk (one call per model); every model in
# RUN_SUMMARY.md went through this same flow.
#
#   1. Register + init (`shortcut gls --unlock`): extract the .rds, register it
#      in the workspace, create the project and write pipeline_config.yaml,
#      WITHOUT computing. Only when no config exists yet: a requeue must not
#      re-init, so the step 1-6 sentinels survive and the pipeline resumes.
#   2. Pin the step-6 / step-7 pools in that config. Serial step 6: the forked
#      ProcessPoolExecutor otherwise deadlocks. Capped step 7: un-pinned,
#      phenoplier sizes the pool from the whole node and thrashes it.
#   3. `workflow gls run` (honours the edited config; resume-safe).
#   4. Gate: refuse to publish a summary that is not COMPLETE
#      (verify_summary.py) -- `summarize` concatenates whatever exists, so an
#      interrupted step 7 yields a smaller but valid-looking file.
#   5. Publish the summary (+ pickle, + the trait filter's exclusion log)
#      atomically next to each other.
#
# Stale projects: phenoplier-cli skips extraction when a model key is already
# registered, and this script skips the run when the project already holds a
# summary -- so a refit .rds would silently reuse the old results. The .rds is
# therefore fingerprinted (sha256) when the project is initialised here, and
# a later call with a different fingerprint removes the registration, moves
# the project aside (<project>.stale.<timestamp>) and starts over. Projects
# that predate this script (no fingerprint on record) are reused as-is, with
# a warning.
#
# phenoplier-cli lives in its own conda env (see setup_env.sh) and `phenoplier`
# shells out to a bare `snakemake`, whose rules shell out to a bare
# `phenoplier`, so the env's bin has to be on PATH.
set -euo pipefail

usage() {
  cat >&2 <<USAGE
usage: $0 --rds <model.rds> --name <project> --summary-out <out.tsv.gz> [options]
  --conda-env NAME            (phenoplier-cli-neo)
  --require-version X.Y.Z     fail unless the installed phenoplier-cli is this release
  --namespace NS              workspace model namespace (clamp)
  --cohort NAME               (phenomexcan_rapid_gwas)
  --lv-percentile FLOAT       (0.01)
  --trait-filter all|biomedical  (biomedical)
  --expected-phenotypes N     completeness gate (2366 = biomedical)
  --executor local|slurm      phenoplier-cli's own executor (local)
  --cluster NAME              phenoplier-cli cluster profile, slurm only ('')
  --n-jobs N                  cores for the local executor (4)
  --workers-step6 N --blas-step6 N --workers-step7 N --blas-step7 N  (1 8 6 2)
Workspace: \$PHENOPLIER_HOME (default \$HOME/phenoplier).
USAGE
  exit 2
}

conda_env=phenoplier-cli-neo
require_version=""
namespace=clamp
cohort=phenomexcan_rapid_gwas
lv_percentile=0.01
trait_filter=biomedical
expected_phenotypes=2366
executor=local
cluster=""
n_jobs=4
workers_step6=1
blas_step6=8
workers_step7=6
blas_step7=2
rds=""
name=""
summary_out=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --rds) rds="$2"; shift 2 ;;
    --name) name="$2"; shift 2 ;;
    --summary-out) summary_out="$2"; shift 2 ;;
    --conda-env) conda_env="$2"; shift 2 ;;
    --require-version) require_version="$2"; shift 2 ;;
    --namespace) namespace="$2"; shift 2 ;;
    --cohort) cohort="$2"; shift 2 ;;
    --lv-percentile) lv_percentile="$2"; shift 2 ;;
    --trait-filter) trait_filter="$2"; shift 2 ;;
    --expected-phenotypes) expected_phenotypes="$2"; shift 2 ;;
    --executor) executor="$2"; shift 2 ;;
    --cluster) cluster="$2"; shift 2 ;;
    --n-jobs) n_jobs="$2"; shift 2 ;;
    --workers-step6) workers_step6="$2"; shift 2 ;;
    --blas-step6) blas_step6="$2"; shift 2 ;;
    --workers-step7) workers_step7="$2"; shift 2 ;;
    --blas-step7) blas_step7="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "unknown argument: $1" >&2; usage ;;
  esac
done
[[ -n "$rds" && -n "$name" && -n "$summary_out" ]] || usage
[[ -f "$rds" ]] || { echo "model not found: $rds" >&2; exit 1; }

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
workspace="${PHENOPLIER_HOME:-$HOME/phenoplier}"
export PHENOPLIER_HOME="$workspace" PHENOPLIER_ROOT_DIR="$workspace"

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$conda_env"
export PATH="${CONDA_PREFIX:-}/bin:${PATH}"

# `shortcut gls --name X` creates the project at exactly <workspace>/projects/X
# (phenoplier-cli: commands/models.py derive_folder_name) -- no date prefix.
project="$workspace/projects/$name"
cfg="$project/pipeline_config.yaml"
results="$project/results/gls/phenoplier"
summary="$results/gls-summary-${cohort}.tsv.gz"
fingerprint_file="$project/clamp_source.sha256"

installed="$(phenoplier --version 2>&1 | tail -1 | sed -E 's/.*v([0-9][0-9A-Za-z.+-]*).*/\1/')"
echo "[start] $(date -Is)  name=$name  rds=$rds  executor=$executor  n_jobs=$n_jobs"
echo "[version] phenoplier-cli $installed  workspace=$workspace"
# Results from different phenoplier-cli releases are not comparable; the
# whole model set is produced on the release phenoplier.yaml pins.
if [[ -n "$require_version" && "$installed" != "$require_version" ]]; then
  echo "[ERROR] phenoplier-cli $installed in env '$conda_env' but $require_version is required" \
    "(workflow/config/phenoplier.yaml: version); rebuild the env with scripts/phenoplier/setup_env.sh" >&2
  exit 1
fi

fingerprint="$(sha256sum "$rds" | cut -d' ' -f1)"
if [[ -d "$project" ]]; then
  if [[ -s "$fingerprint_file" ]]; then
    recorded="$(cut -d' ' -f1 < "$fingerprint_file")"
    if [[ "$recorded" != "$fingerprint" ]]; then
      stale="${project}.stale.$(date +%Y%m%dT%H%M%S)"
      echo "[stale] $(date -Is)  model changed since this project was initialised:"
      echo "[stale]   recorded $recorded"
      echo "[stale]   current  $fingerprint"
      echo "[stale]   moving $project -> $stale and re-registering $namespace/$name"
      exec 9>"$workspace/.clamp_registry.lock"; flock 9
      phenoplier model remove --force "$namespace/$name" || true
      flock -u 9; exec 9>&-
      mv "$project" "$stale"
    fi
  elif [[ -s "$cfg" ]]; then
    echo "[warn] $name: project predates fingerprinting; reusing it as-is." \
      "Delete $project to recompute from the current model." >&2
  fi
fi

if [[ ! -s "$summary" ]]; then
  if [[ ! -s "$cfg" ]]; then
    cluster_args=()
    if [[ -n "$cluster" ]]; then
      cluster_args=(--cluster "$cluster")
    fi
    # A fresh project always gets a fresh registration: phenoplier-cli skips
    # extraction when the key already exists, which would keep an old Z.
    # The workspace registry (models/registry.toml) is rewritten whole by
    # `model remove` without a lock, so concurrent GLS jobs initialising at
    # the same time could drop each other's entries: hold a workspace-wide
    # lock across remove + register/init. Only this short phase is
    # serialised; the hours-long pipeline run is not.
    lock="$workspace/.clamp_registry.lock"
    exec 9>"$lock"
    echo "[init] $(date -Is)  waiting for registry lock $lock"
    flock 9
    phenoplier model remove --force "$namespace/$name" 2>/dev/null || true
    echo "[init] $(date -Is)  registering model + writing $cfg"
    phenoplier shortcut gls \
      --input "$rds" \
      --name "$name" \
      --model-namespace "$namespace" \
      --cohort "$cohort" \
      --lv-percentile "$lv_percentile" \
      --trait-filter "$trait_filter" \
      --executor "$executor" \
      --n-jobs "$n_jobs" \
      ${cluster_args[@]+"${cluster_args[@]}"} \
      --unlock
    printf '%s  %s\n' "$fingerprint" "$rds" > "$fingerprint_file"
    flock -u 9
    exec 9>&-
  else
    echo "[init] $(date -Is)  resuming: $cfg exists"
  fi

  python - "$cfg" "$workers_step6" "$blas_step6" "$workers_step7" "$blas_step7" <<'PY'
import sys, yaml
path, nw6, bl6, nw7, bl7 = sys.argv[1], *map(int, sys.argv[2:6])
with open(path) as fh:
    cfg = yaml.safe_load(fh)
cfg["n_workers_step6"] = nw6
cfg["blas_threads_step6"] = bl6
cfg["n_workers_step7"] = nw7
cfg["blas_threads_step7"] = bl7
with open(path, "w") as fh:
    yaml.safe_dump(cfg, fh, sort_keys=False)
print(f"[pin] step6={nw6}x{bl6} step7={nw7}x{bl7} -> {path}")
PY

  rm -f "$project/.snakemake/locks/"*.lock 2>/dev/null || true
  echo "[run] $(date -Is)"
  phenoplier workflow gls run --config-file "$cfg"
else
  echo "[run] skip: $summary exists"
fi

# The summary is named after the cohort (gls-summary-<cohort-slug>.tsv.gz);
# older phenoplier-cli kept a gls-summary-phenomexcan.* alias. Take the
# cohort-named file, then the alias, then whatever single summary exists.
if [[ ! -s "$summary" ]]; then
  summary="$results/gls-summary-phenomexcan.tsv.gz"
fi
if [[ ! -s "$summary" ]]; then
  summary="$(ls -1 "$results"/gls-summary-*.tsv.gz 2>/dev/null | head -n1 || true)"
fi
if [[ -z "$summary" || ! -s "$summary" ]]; then
  echo "[ERROR] $name: no gls-summary-*.tsv.gz under $results (run failed?)" >&2
  exit 1
fi

# Gate: complete (expected phenotype count), no zero p-values, NaN only on
# LVs flagged degenerate. A failed gate publishes nothing AND drops the
# project's combined summary, so the next attempt re-enters the pipeline
# (resume-safe: step 7 finishes what is missing, summarize reruns) instead of
# skipping straight back to the same incomplete file.
if ! python "$here/verify_summary.py" "$summary" "$expected_phenotypes"; then
  echo "[ERROR] $name: summary failed the completeness gate; removing it so a requeue resumes step 7" >&2
  rm -f "$summary" "${summary%.tsv.gz}.pkl.gz"
  exit 1
fi

mkdir -p "$(dirname "$summary_out")"
tmp="${summary_out}.tmp.$$"
cp -f "$summary" "$tmp"
mv -f "$tmp" "$summary_out"
stem="${summary_out%.tsv.gz}"
if [[ -f "${summary%.tsv.gz}.pkl.gz" ]]; then
  cp -f "${summary%.tsv.gz}.pkl.gz" "${stem}.pkl.gz"
fi
# The filter's exclusion log travels with the results: which traits were
# dropped, and why, is part of the provenance of these numbers.
if [[ -f "$project/trait_filter_excluded.tsv" ]]; then
  cp -f "$project/trait_filter_excluded.tsv" "${stem}_trait_filter_excluded.tsv"
fi
echo "[done] $(date -Is)  $name -> $summary_out"
