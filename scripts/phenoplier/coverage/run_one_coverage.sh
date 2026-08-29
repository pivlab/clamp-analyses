#!/usr/bin/env bash
# Run `phenoplier shortcut gls` + `store build` for ONE ARCHS4 coverage model
# (fraction rs, seed). Idempotent: skips any stage whose output already exists,
# so the worklist can be safely re-launched. Streams the model from the source
# host if not already staged locally.
#
# Usage: run_one_coverage.sh <rs> <seed>
# All paths env-overridable (defaults = the local /media/data setup).
set -euo pipefail

rs="$1"; seed="$2"

ENVBIN="${PHENOPLIER_ENVBIN:-/home/haoyu/miniforge3/envs/phenoplier-cli-neo/bin}"
WS="${PHENOPLIER_HOME:-/media/data/phenoplier_workspace}"
BASEDIR="${CLAMP_COV_DIR:-/media/data/clamp_coverage}"
SRC_SSH="${CLAMP_COV_SRC_SSH:-pico}"
SRC_ROOT="${CLAMP_COV_SRC_ROOT:-/pividori_lab/marc_projects/clamp-analyses-coverage-bp/output/01_model_building/02_archs4/02_coverage_study/models/archs4}"
COHORT="${CLAMP_COHORT:-phenomexcan_rapid_gwas}"
TRAIT_FILTER="${CLAMP_TRAIT_FILTER:-biomedical}"
NJOBS="${CLAMP_NJOBS:-24}"

export PATH="$ENVBIN:$PATH"
export PHENOPLIER_HOME="$WS"
export PHENOPLIER_ROOT_DIR="$WS"
# Cap per-process BLAS/OpenMP threads to 1 so parallelism comes from the JOB
# level (step 3 runs 22 chromosomes; step 6/7 run worker pools) rather than each
# job also spawning a full set of BLAS threads. Without this, on a workstation
# with no SLURM cgroup a single model oversubscribes the box (22 chroms x N BLAS
# threads -> load >80 on 24 cores) and crawls. SLURM sites get this via -c.
export OMP_NUM_THREADS="${CLAMP_BLAS_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${CLAMP_BLAS_THREADS:-1}"
export MKL_NUM_THREADS="${CLAMP_BLAS_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${CLAMP_BLAS_THREADS:-1}"
export VECLIB_MAXIMUM_THREADS="${CLAMP_BLAS_THREADS:-1}"

name="cov_rs${rs}_seed${seed}_CLAMPfull_bp"
stage="$BASEDIR/models/rs${rs}_seed${seed}"
rds="$stage/CLAMPfull_bp.rds"
log="$BASEDIR/logs/${name}.log"
summary="$WS/projects/$name/results/gls/phenoplier/gls-summary-phenomexcan.tsv.gz"
store="$BASEDIR/stores/${name}.h5"
mkdir -p "$stage" "$BASEDIR/stores" "$BASEDIR/logs"

{
  echo "===== $name  $(date -Is)  (affinity: $(taskset -pc $$ 2>/dev/null | sed 's/.*: //')) ====="
  if [ ! -s "$rds" ]; then
    echo "[stream] ${SRC_SSH}:${SRC_ROOT}/rs${rs}/seed${seed} -> $rds  $(date -Is)"
    ssh "$SRC_SSH" "cat ${SRC_ROOT}/rs${rs}/seed${seed}/CLAMPfull_bp/CLAMPfull_bp.rds" > "$rds"
  fi
  if [ ! -s "$summary" ]; then
    echo "[gls] $(date -Is)  --trait-filter ${TRAIT_FILTER} --n-jobs ${NJOBS}"
    phenoplier shortcut gls --input "$rds" --name "$name" --model-namespace clamp \
      --trait-filter "$TRAIT_FILTER" --lv-percentile 0.01 --executor local --n-jobs "$NJOBS"
  else
    echo "[gls] skip (summary exists)"
  fi
  if [ ! -s "$store" ]; then
    echo "[store] $(date -Is)"
    phenoplier store build --output "$store" --cohort "$COHORT" --clamp-rds "$rds" \
      --gls-dir "$WS/projects/$name/results/gls/phenoplier/${COHORT}" --method-name CLAMP
  else
    echo "[store] skip (exists)"
  fi
  echo "[done] $name  $(date -Is)"
} >> "$log" 2>&1
