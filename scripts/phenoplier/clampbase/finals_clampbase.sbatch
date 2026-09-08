#!/bin/bash
#SBATCH -J cb_final
#SBATCH -p defq
#SBATCH -c 12
#SBATCH -o /pividori_lab/phenoplier_workspace/clampbase_finals/logs/%x_%A_%a.out
#
# Full CLAMPbase finals (archs4/gtex/recount2). Submit per memory tier:
#   sbatch --mem=88G --array=0   finals_clampbase.sbatch   # archs4 (7.8G rds)
#   sbatch --mem=48G --array=1-2 finals_clampbase.sbatch   # gtex, recount2
set -uo pipefail
TISSUE=(archs4 gtex recount2)
i="${SLURM_ARRAY_TASK_ID:?run as --array=0-2}"
t="${TISSUE[$i]}"

source ~/miniforge3/etc/profile.d/conda.sh
conda activate phenoplier-cli-neo
export PATH="$CONDA_PREFIX/bin:$PATH"
export PHENOPLIER_HOME=/pividori_lab/phenoplier_workspace
export CLAMPBASE_CENTRAL=/pividori_lab/phenoplier_workspace/clampbase_finals
mkdir -p "$CLAMPBASE_CENTRAL/logs"

rds="/pividori_lab/marc_projects/clamp-analyses/output/98_final_models/clampbase/${t}/CLAMPbase.rds"
name="final_${t}_CLAMPbase"
bash "${CLAMPBASE_RUNONE:-/pividori_lab/haoyu_projects/clampbase_scripts/_run_clampbase_model.sh}" "$name" "$rds"
