#!/bin/bash
#SBATCH -J cov_gls
#SBATCH -p defq
#SBATCH -c 32
#SBATCH --mem=150G
#SBATCH -o /pividori_lab/phenoplier_workspace/clamp_coverage/logs/%x_%A_%a.out
#
# Run the BIG coverage models on pico (defq): rs50 s2/s3, rs75 s1/s2/s3,
# rs100 s2/s3 -- the 7 that were queued on the local box. pico is faster and its
# saturation run is finishing, so these backfill around it. Submit throttled:
#   sbatch --array=0-6%2 pico_coverage.sbatch
# SLURM won't oversubscribe the node, so %2 (2 x 32 = 64 cores) runs alongside
# the finishing saturation jobs and grows into freed cores. Idempotent +
# resume-safe. The coverage models live on pico's filesystem -- no transfer.
set -euo pipefail
RS=(50 50 75 75 75 100 100); SEED=(2 3 1 2 3 2 3)
i="${SLURM_ARRAY_TASK_ID:?run as --array=0-6%2}"
rs="${RS[$i]}"; seed="${SEED[$i]}"

WS=/pividori_lab/phenoplier_workspace
SRC=/pividori_lab/marc_projects/clamp-analyses-coverage-bp/output/01_model_building/02_archs4/02_coverage_study/models/archs4
COHORT=phenomexcan_rapid_gwas

source ~/miniforge3/etc/profile.d/conda.sh
conda activate phenoplier-cli-neo
export PATH="$CONDA_PREFIX/bin:$PATH"
export PHENOPLIER_HOME="$WS" PHENOPLIER_ROOT_DIR="$WS"

name="cov_rs${rs}_seed${seed}_CLAMPfull_bp"
rds="$SRC/rs${rs}/seed${seed}/CLAMPfull_bp/CLAMPfull_bp.rds"
summary="$WS/projects/$name/results/gls/phenoplier/gls-summary-phenomexcan_rapid_gwas.tsv.gz"
store="$WS/clamp_coverage/stores/${name}.h5"
NJOBS="${SLURM_CPUS_PER_TASK:-32}"
mkdir -p "$WS/clamp_coverage/stores" "$WS/clamp_coverage/logs"

echo "[start] $(date -Is)  $name  njobs=$NJOBS  ver=$(phenoplier -v 2>&1 | tail -1)"
[ -s "$summary" ] || rm -f "$WS/projects/$name/.snakemake/locks/"*.lock 2>/dev/null || true
[ -s "$summary" ] || phenoplier shortcut gls --input "$rds" --name "$name" \
  --model-namespace clamp --trait-filter biomedical --lv-percentile 0.01 \
  --executor local --n-jobs "$NJOBS"
[ -s "$store" ] || phenoplier store build --output "$store" --cohort "$COHORT" \
  --clamp-rds "$rds" --gls-dir "$WS/projects/$name/results/gls/phenoplier/${COHORT}" \
  --method-name CLAMP
echo "[done] $(date -Is)  $name"
