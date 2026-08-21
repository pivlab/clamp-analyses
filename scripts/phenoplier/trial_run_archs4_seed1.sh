#!/usr/bin/env bash
# Trial end-to-end run of the phenoplier-cli integration against a real
# production CLAMP model (ARCHS4 final, seed 1, 1,728 LVs), documented in
# nbs/03_model_biology/00_phenoplier/00_trial_run_archs4_final_seed1.ipynb.
#
# This is deliberately a direct `phenoplier shortcut gls` call rather than a
# Snakemake rule: it validates the CLI itself, before the DAG in
# workflow/rules/phenoplier.smk drives it across the whole model matrix.
#
#   bash scripts/phenoplier/trial_run_archs4_seed1.sh
#
# Long-running (steps 1-5 take many hours regardless of LV count); run it
# detached and follow the log:
#
#   setsid nohup bash scripts/phenoplier/trial_run_archs4_seed1.sh \
#       > trial_run.log 2>&1 < /dev/null &
set -euo pipefail

CONDA_ENV="${PHENOPLIER_ENV:-phenoplier-cli-neo}"
export PHENOPLIER_ROOT_DIR="${PHENOPLIER_ROOT_DIR:-/pividori_lab/phenoplier_workspace}"

MODEL="${MODEL_RDS:-/pividori_lab/marc_projects/clamp-analyses-coverage-bp/output/01_model_building/02_archs4/01_final_model/seed1/CLAMPfull.rds}"
NAME="${MODEL_NAME:-archs4_final_seed1_CLAMPfull}"
TRAIT_FILTER="${TRAIT_FILTER:-biomedical}"
N_JOBS="${N_JOBS:-22}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

# `phenoplier` shells out to a bare `snakemake`, and every Snakemake rule
# shells out to a bare `phenoplier`: the env's bin must be ON PATH. Calling
# the executables by absolute path is not enough.
export PATH="${CONDA_PREFIX}/bin:${PATH}"

echo "[start]     $(date -Is)"
echo "[version]   $(phenoplier -v 2>&1 | tail -1)"
echo "[model]     ${MODEL}"
echo "[filter]    ${TRAIT_FILTER}"
echo "[workspace] ${PHENOPLIER_ROOT_DIR}"

time phenoplier shortcut gls \
  --input "${MODEL}" \
  --name "${NAME}" \
  --model-namespace clamp \
  --lv-percentile 0.01 \
  --trait-filter "${TRAIT_FILTER}" \
  --executor local \
  --n-jobs "${N_JOBS}"

echo "[done]      $(date -Is)"
