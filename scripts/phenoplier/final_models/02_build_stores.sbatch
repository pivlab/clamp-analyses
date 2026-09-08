#!/bin/bash
#SBATCH -J clamp_store
#SBATCH -c 4
#SBATCH --mem=96G
#SBATCH -o %x_%j.out
#
# Step 3 of the final-models GLS run (PR #28, item 7). Composite each final
# model's CLAMP .rds (Z matrix, hyperparameters, provenance) with that model's
# GLS phenomexcan results into a one-file HDF5 composite study store, via
# `phenoplier store build` (phenoplier-cli #70). One store per model, because a
# store holds a single /clamp model but may hold many GWAS cohorts.
# Run after the GLS jobs finish:  sbatch 02_build_stores.sbatch
#
# Requires `anndata` (+ zarr) in the env -- the store's element encodings use
# anndata.io. On an air-gapped cluster (pico has no PyPI), transfer the wheels
# from a networked host and `pip install --no-index --no-deps` them; keep the
# pinned numpy/pandas untouched. See README.
set -euo pipefail

CONDA_SH="${PHENOPLIER_CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
ENV_NAME="${PHENOPLIER_ENV:-phenoplier-cli-neo}"
WS="${PHENOPLIER_HOME:-/pividori_lab/phenoplier_workspace}"
MODELS_DIR="${CLAMP_FINAL_MODELS:-/pividori_lab/marc_projects/clamp-analyses-coverage-bp/output/98_final_models}"
OUT_DIR="${CLAMP_STORE_DIR:-$WS/clamp_final_stores}"
DATASETS="${CLAMP_FINAL_DATASETS:-recount2 gtex archs4}"
COHORT="${CLAMP_COHORT:-phenomexcan_rapid_gwas}"

# shellcheck disable=SC1091
source "$CONDA_SH"
conda activate "$ENV_NAME"
export PATH="$CONDA_PREFIX/bin:$PATH"   # store build --clamp-rds needs R on PATH
export PHENOPLIER_HOME="$WS"
export PHENOPLIER_ROOT_DIR="$WS"
mkdir -p "$OUT_DIR"

echo "[start] $(date -Is)  host=$(hostname)  version=$(phenoplier -v 2>&1 | tail -1)"
for DS in $DATASETS; do
  proj="$WS/projects/final_${DS}_CLAMPfull_bp"
  out="$OUT_DIR/final_${DS}_CLAMPfull_bp.h5"
  echo "===== store $DS -> $out ====="
  phenoplier store build \
    --output "$out" \
    --cohort "$COHORT" \
    --clamp-rds "$MODELS_DIR/${DS}/CLAMPfull_bp.rds" \
    --gls-dir "$proj/results/gls/phenoplier/${COHORT}" \
    --method-name CLAMP
done
echo "[done] $(date -Is)"
