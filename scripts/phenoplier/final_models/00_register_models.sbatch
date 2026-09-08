#!/bin/bash
#SBATCH -J clamp_reg
#SBATCH -c 8
#SBATCH --mem=160G
#SBATCH -o %x_%j.out
#
# Step 1 of the final-models GLS run (PR #28, item 0). Register the three CLAMP
# final models (one CLAMPfull_bp per dataset) into the workspace ONCE,
# sequentially. The GLS jobs (01_run_final_gls.sbatch) then reference them with
# `shortcut gls --model`, which only READS the shared models/registry.toml -- so
# the three GLS jobs can run in parallel without racing on registration writes
# (registry.toml is not written atomically). Gate them with
#   sbatch --dependency=afterok:<this_jobid> 01_run_final_gls.sbatch <dataset>
#
# Partition/account/QOS/time-limit come from `sbatch` CLI flags or cluster
# defaults, not baked in (see README for pico vs Alpine). Paths default to the
# pico layout where this was run; override via the env vars below.
set -euo pipefail

CONDA_SH="${PHENOPLIER_CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
ENV_NAME="${PHENOPLIER_ENV:-phenoplier-cli-neo}"
WS="${PHENOPLIER_HOME:-/pividori_lab/phenoplier_workspace}"
MODELS_DIR="${CLAMP_FINAL_MODELS:-/pividori_lab/marc_projects/clamp-analyses-coverage-bp/output/98_final_models}"
DATASETS="${CLAMP_FINAL_DATASETS:-archs4 gtex recount2}"
NAMESPACE="${PHENOPLIER_NAMESPACE:-clamp}"

# shellcheck disable=SC1091
source "$CONDA_SH"
conda activate "$ENV_NAME"
export PATH="$CONDA_PREFIX/bin:$PATH"
export PHENOPLIER_HOME="$WS"
export PHENOPLIER_ROOT_DIR="$WS"

echo "[start] $(date -Is)  host=$(hostname)  version=$(phenoplier -v 2>&1 | tail -1)"
for DS in $DATASETS; do
  echo "=== register $DS ==="
  phenoplier model register "$MODELS_DIR/${DS}/CLAMPfull_bp.rds" \
    --name "final_${DS}_CLAMPfull_bp" --namespace "$NAMESPACE"
done
echo "=== registered ==="
phenoplier model list 2>&1 | grep -i "final_.*CLAMPfull_bp" || true
echo "[done] $(date -Is)"
