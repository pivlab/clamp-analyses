#!/bin/bash
#SBATCH -J cb_cov
#SBATCH -p defq
#SBATCH -c 12
#SBATCH -o /pividori_lab/phenoplier_workspace/clampbase_coverage/logs/%x_%A_%a.out
#
# Coverage-study CLAMPbase (21 models). Memory tiered by fraction; submit per
# tier (indices are grouped by rs, smallest first):
#   sbatch --mem=48G --array=0-8%14   cov_clampbase.sbatch   # rs1, rs5, rs10
#   sbatch --mem=56G --array=9-14%14  cov_clampbase.sbatch   # rs25, rs50
#   sbatch --mem=64G --array=15-17%14 cov_clampbase.sbatch   # rs75
#   sbatch --mem=88G --array=18-20%14 cov_clampbase.sbatch   # rs100 (full)
set -uo pipefail
RS=(1 1 1  5 5 5  10 10 10  25 25 25  50 50 50  75 75 75  100 100 100)
SD=(1 2 3  1 2 3  1  2  3   1  2  3   1  2  3   1  2  3   1   2   3)
i="${SLURM_ARRAY_TASK_ID:?run as --array=0-20}"
rs="${RS[$i]}"; seed="${SD[$i]}"

# fraction -> level dir
declare -A LVL=( [1]=00_bp_coverage_study_01 [5]=01_bp_coverage_study_05 \
  [10]=02_bp_coverage_study_10 [25]=03_bp_coverage_study_25 \
  [50]=04_bp_coverage_study_50 [75]=05_bp_coverage_study_75 \
  [100]=06_bp_coverage_study_100 )

source ~/miniforge3/etc/profile.d/conda.sh
conda activate phenoplier-cli-neo
export PATH="$CONDA_PREFIX/bin:$PATH"
export PHENOPLIER_HOME=/pividori_lab/phenoplier_workspace
export CLAMPBASE_CENTRAL=/pividori_lab/phenoplier_workspace/clampbase_coverage
mkdir -p "$CLAMPBASE_CENTRAL/logs"

COV=/pividori_lab/marc_projects/clamp-analyses-coverage-bp/output/01_model_building/02_archs4/02_coverage_study
rds="$COV/${LVL[$rs]}/study_coverage_rs${rs}_seed_${seed}/CLAMPbase.rds"
name="cov_rs${rs}_seed${seed}_CLAMPbase"
bash "${CLAMPBASE_RUNONE:-/pividori_lab/haoyu_projects/clampbase_scripts/_run_clampbase_model.sh}" "$name" "$rds"
