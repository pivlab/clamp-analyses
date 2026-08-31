#!/bin/bash
#SBATCH -J cov_gls
#SBATCH -p acpu
#SBATCH --account=amc-general
#SBATCH --qos=cpu-normal
#SBATCH -c 16
#SBATCH --mem=80G
#SBATCH --time=16:00:00
#SBATCH -o /pl/active/pivlab/projects/hzhang/clamp_coverage/logs/%x_%A_%a.out
#
# Alpine backfill share of the ARCHS4 coverage study: the 6 smallest models
# (rs1×3, rs5×3), one per array task. GLS + store, idempotent. Submit:
#   sbatch --array=0-5%6 alpine_coverage.sbatch
# (SLURM cgroup caps step-6/7 pools via -c, so no taskset needed here.)
set -euo pipefail
RS=(1 1 1 5 5 5); SEED=(1 2 3 1 2 3)
i="${SLURM_ARRAY_TASK_ID:?run as --array=0-5%6}"
rs="${RS[$i]}"; seed="${SEED[$i]}"

MF=/pl/active/pivlab/projects/hzhang/miniforge3
PREFIX=/projects/hzhang24@xsede.org/software/anaconda/envs/phenoplier-cli-neo
WS=/pl/active/pivlab/projects/hzhang/phenoplier_workspace
BASE=/pl/active/pivlab/projects/hzhang/clamp_coverage
COHORT=phenomexcan_rapid_gwas

source "$MF/etc/profile.d/conda.sh"
conda activate "$PREFIX"
export PATH="$PREFIX/bin:$PATH"
export PHENOPLIER_HOME="$WS" PHENOPLIER_ROOT_DIR="$WS"

name="cov_rs${rs}_seed${seed}_CLAMPfull_bp"
rds="$BASE/models/rs${rs}_seed${seed}/CLAMPfull_bp.rds"
summary="$WS/projects/$name/results/gls/phenoplier/gls-summary-phenomexcan.tsv.gz"
store="$BASE/stores/${name}.h5"
NJOBS="${SLURM_CPUS_PER_TASK:-16}"
mkdir -p "$BASE/stores"

echo "[start] $(date -Is)  $name  njobs=$NJOBS  ver=$(phenoplier -v 2>&1 | tail -1)"
# Resume-safe: clear a stale Snakemake lock from a prior timed-out run so a
# re-submitted task resumes from its last completed step, not from scratch.
[ -s "$summary" ] || rm -f "$WS/projects/$name/.snakemake/locks/"*.lock 2>/dev/null || true
[ -s "$summary" ] || phenoplier shortcut gls --input "$rds" --name "$name" \
  --model-namespace clamp --trait-filter biomedical --lv-percentile 0.01 \
  --executor local --n-jobs "$NJOBS"
[ -s "$store" ] || phenoplier store build --output "$store" --cohort "$COHORT" \
  --clamp-rds "$rds" --gls-dir "$WS/projects/$name/results/gls/phenoplier/${COHORT}" \
  --method-name CLAMP
echo "[done] $(date -Is)  $name"
