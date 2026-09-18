#!/bin/bash
#SBATCH -J sat_reg
#SBATCH -c 8
#SBATCH --mem=160G
#SBATCH -o %x_%j.out
#
# Step 1 of the ARCHS4 saturation k=1728 GLS run (stacked on clamp-analyses
# PR #28). Register the 17 saturation models ONCE, sequentially, so the GLS
# array can reference them with `shortcut gls --model` (read-only on the shared
# registry -- no race, and no re-extracting each up-to-5.7 GB .rds per job).
# Gate the GLS array on this with --dependency=afterok:<this_jobid>.
#
# Models: rs{1,5,10,25,50,75} x k1728 x seed{1,2,3}; rs1/seed1 is missing on
# pico, so 17 not 18. Paths default to pico; override via env vars.
set -euo pipefail

CONDA_SH="${PHENOPLIER_CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
ENV_NAME="${PHENOPLIER_ENV:-phenoplier-cli-neo}"
WS="${PHENOPLIER_HOME:-/pividori_lab/phenoplier_workspace}"
SAT_ROOT="${CLAMP_SAT_ROOT:-/pividori_lab/marc_projects/clamp-analyses-coverage-bp/output/01_model_building/02_archs4/03_saturation_study/models}"
NAMESPACE="${PHENOPLIER_NAMESPACE:-clamp}"

# 17 (rs, seed) cells, big models first (rs75 -> rs1).
RS=(75 75 75 50 50 50 25 25 25 10 10 10 5 5 5 1 1)
SEED=(1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 2 3)

# shellcheck disable=SC1091
source "$CONDA_SH"
conda activate "$ENV_NAME"
export PATH="$CONDA_PREFIX/bin:$PATH"
export PHENOPLIER_HOME="$WS"
export PHENOPLIER_ROOT_DIR="$WS"

echo "[start] $(date -Is)  host=$(hostname)  version=$(phenoplier -v 2>&1 | tail -1)"
for i in "${!RS[@]}"; do
  rs="${RS[$i]}"; seed="${SEED[$i]}"
  rds="$SAT_ROOT/rs${rs}/k1728/seed${seed}/CLAMPfull_bp/CLAMPfull_bp.rds"
  name="sat_rs${rs}_k1728_seed${seed}_CLAMPfull_bp"
  echo "=== register [$i] $name ==="
  phenoplier model register "$rds" --name "$name" --namespace "$NAMESPACE"
done
echo "=== registered ==="
phenoplier model list 2>&1 | grep -i "sat_rs.*k1728" || true
echo "[done] $(date -Is)"
