#!/bin/bash
#SBATCH -J cov_tiny
#SBATCH -p defq
#SBATCH -c 22
#SBATCH --mem=100G
#SBATCH -o /pividori_lab/phenoplier_workspace/clamp_coverage/logs/%x_%A_%a.out
#
# The 3 tiny coverage models that kept stalling on Alpine (ProcessPoolExecutor
# deadlock at step 6): rs1 s1/s3, rs5 s2. Run them fresh on pico (fast, reliable,
# models on-disk). Submit:  sbatch --array=0-2%3 pico_coverage_tiny.sbatch
set -euo pipefail
RS=(1 1 5); SEED=(1 3 2)
i="${SLURM_ARRAY_TASK_ID:?run as --array=0-2%3}"
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
NJOBS="${SLURM_CPUS_PER_TASK:-22}"
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
