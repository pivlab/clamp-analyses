#!/bin/bash
#SBATCH -J sat_store
#SBATCH -c 4
#SBATCH --mem=64G
#SBATCH -o %x_%A_%a.out
#
# Step 3: build one HDF5 composite study store per saturation model (CLAMP .rds
# + that model's GLS results), as a %6 array after the GLS array finishes:
#   sbatch --dependency=afterok:<gls_array_jobid> --array=0-16%6 02_build_stores.sbatch
set -euo pipefail

CONDA_SH="${PHENOPLIER_CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
ENV_NAME="${PHENOPLIER_ENV:-phenoplier-cli-neo}"
WS="${PHENOPLIER_HOME:-/pividori_lab/phenoplier_workspace}"
SAT_ROOT="${CLAMP_SAT_ROOT:-/pividori_lab/marc_projects/clamp-analyses-coverage-bp/output/01_model_building/02_archs4/03_saturation_study/models}"
OUT_DIR="${CLAMP_SAT_STORE_DIR:-$WS/clamp_saturation_stores}"
COHORT="${CLAMP_COHORT:-phenomexcan_rapid_gwas}"

RS=(75 75 75 50 50 50 25 25 25 10 10 10 5 5 5 1 1)
SEED=(1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 2 3)

i="${SLURM_ARRAY_TASK_ID:?run as a job array (--array=0-16%6)}"
rs="${RS[$i]}"; seed="${SEED[$i]}"
name="sat_rs${rs}_k1728_seed${seed}_CLAMPfull_bp"

# shellcheck disable=SC1091
source "$CONDA_SH"
conda activate "$ENV_NAME"
export PATH="$CONDA_PREFIX/bin:$PATH"   # store build --clamp-rds needs R on PATH
export PHENOPLIER_HOME="$WS"
export PHENOPLIER_ROOT_DIR="$WS"
mkdir -p "$OUT_DIR"

rds="$SAT_ROOT/rs${rs}/k1728/seed${seed}/CLAMPfull_bp/CLAMPfull_bp.rds"
proj="$WS/projects/${name}"
out="$OUT_DIR/${name}.h5"
echo "[start] $(date -Is)  task=$i  $name -> $out"
phenoplier store build \
  --output "$out" \
  --cohort "$COHORT" \
  --clamp-rds "$rds" \
  --gls-dir "$proj/results/gls/phenoplier/${COHORT}" \
  --method-name CLAMP
echo "[done] $(date -Is)  task=$i"
