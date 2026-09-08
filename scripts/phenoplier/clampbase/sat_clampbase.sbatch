#!/bin/bash
#SBATCH -J cb_sat
#SBATCH -p defq
#SBATCH -c 12
#SBATCH --mem=64G
#SBATCH -o /pividori_lab/phenoplier_workspace/clampbase_saturation/logs/%x_%A_%a.out
#
# Saturation-study CLAMPbase, k1728 only (17 models; rs1_seed1 has no model).
# Uniform 64G. Submit:  sbatch --array=0-16%14 sat_clampbase.sbatch
set -uo pipefail
# 17 (rs,seed) pairs, rs1_seed1 excluded.
RS=(1 1  5 5 5  10 10 10  25 25 25  50 50 50  75 75 75)
SD=(2 3  1 2 3  1  2  3   1  2  3   1  2  3   1  2  3)
i="${SLURM_ARRAY_TASK_ID:?run as --array=0-16}"
rs="${RS[$i]}"; seed="${SD[$i]}"

source ~/miniforge3/etc/profile.d/conda.sh
conda activate phenoplier-cli-neo
export PATH="$CONDA_PREFIX/bin:$PATH"
export PHENOPLIER_HOME=/pividori_lab/phenoplier_workspace
export CLAMPBASE_CENTRAL=/pividori_lab/phenoplier_workspace/clampbase_saturation
mkdir -p "$CLAMPBASE_CENTRAL/logs"

SAT=/pividori_lab/marc_projects/clamp-analyses-coverage-bp/output/01_model_building/02_archs4/03_saturation_study
rds="$SAT/study_saturation_rs${rs}_k1728_seed_${seed}/CLAMPbase.rds"
name="sat_rs${rs}_k1728_seed${seed}_CLAMPbase"
bash "${CLAMPBASE_RUNONE:-/pividori_lab/haoyu_projects/clampbase_scripts/_run_clampbase_model.sh}" "$name" "$rds"
