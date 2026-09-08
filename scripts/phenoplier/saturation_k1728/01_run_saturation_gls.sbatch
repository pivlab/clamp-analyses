#!/bin/bash
#SBATCH -J sat_gls
#SBATCH -c 22
#SBATCH --mem=100G
#SBATCH -o %x_%A_%a.out
#
# Step 2: run `phenoplier shortcut gls` for the 17 saturation k=1728 models as a
# SLURM job array THROTTLED to 6 concurrent, so we stay under 80% of pico
# (6 x 22 = 132 of 168 cores = 78.6%, 6 x 100G = 600G of 900G). Submit:
#   REG=$(sbatch --parsable 00_register_models.sbatch)
#   sbatch --dependency=afterok:$REG --array=0-16%6 01_run_saturation_gls.sbatch
#
# Array indices are big-model-first (rs75 -> rs1), so the three long rs75 jobs
# start in the first wave and the smaller ones backfill. `-c 22` covers step 3's
# 22 chromosomes in one wave; phenoplier sizes its step-6/7 pools from the
# cgroup (sched_getaffinity), so this both caps and pins them.
set -euo pipefail

CONDA_SH="${PHENOPLIER_CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
ENV_NAME="${PHENOPLIER_ENV:-phenoplier-cli-neo}"
WS="${PHENOPLIER_HOME:-/pividori_lab/phenoplier_workspace}"
NAMESPACE="${PHENOPLIER_NAMESPACE:-clamp}"
TRAIT_FILTER="${CLAMP_TRAIT_FILTER:-biomedical}"
LV_PERCENTILE="${CLAMP_LV_PERCENTILE:-0.01}"

RS=(75 75 75 50 50 50 25 25 25 10 10 10 5 5 5 1 1)
SEED=(1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 2 3)

i="${SLURM_ARRAY_TASK_ID:?run as a job array (--array=0-16%6)}"
rs="${RS[$i]}"; seed="${SEED[$i]}"
model="${NAMESPACE}/sat_rs${rs}_k1728_seed${seed}_CLAMPfull_bp"

# shellcheck disable=SC1091
source "$CONDA_SH"
conda activate "$ENV_NAME"
export PATH="$CONDA_PREFIX/bin:$PATH"
export PHENOPLIER_HOME="$WS"
export PHENOPLIER_ROOT_DIR="$WS"

NJOBS="${SLURM_CPUS_PER_TASK:-16}"
echo "[start] $(date -Is)  host=$(hostname)  task=$i  model=$model  njobs=$NJOBS  filter=$TRAIT_FILTER"
echo "[version] $(phenoplier -v 2>&1 | tail -1)"

time phenoplier shortcut gls \
  --model "$model" \
  --lv-percentile "$LV_PERCENTILE" \
  --trait-filter "$TRAIT_FILTER" \
  --executor local \
  --n-jobs "$NJOBS"
echo "[done] $(date -Is)  task=$i"
