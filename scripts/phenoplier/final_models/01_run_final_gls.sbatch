#!/bin/bash
#SBATCH -J clamp_gls
#SBATCH -c 32
#SBATCH --mem=150G
#SBATCH -o %x_%j.out
#
# Step 2 of the final-models GLS run (PR #28, item 0). Run `phenoplier shortcut
# gls` for ONE pre-registered CLAMP final model, on a single node with the
# LOCAL executor. Submit one job per dataset; they run in parallel:
#   R=$(sbatch --parsable 00_register_models.sbatch)
#   for ds in archs4 gtex recount2; do
#     sbatch --dependency=afterok:$R 01_run_final_gls.sbatch $ds
#   done
#
# `-c 32` gives the step-6/7 process pools a 32-CPU cgroup; phenoplier sizes
# them from sched_getaffinity, so this both caps and pins them correctly (an
# un-cgrouped run auto-detects the whole node and can exhaust it -- the origin
# of the step-7 BrokenProcessPool seen on pico's 168-core node).
#
# Why local-executor rather than `--executor slurm --cluster ...`: `shortcut
# gls` (phenoplier-cli main @ 3394d44) has no --cluster/--partition/--qos
# passthrough to `workflow gls init`, so its slurm path can't pick up a site
# profile. Running the whole pipeline locally in one allocation matches the
# validated path. Steps 1-5 (SNP covariance, GWAS/TWAS compile, per-chromosome
# gene correlations) scale with the model's GENE count, not its LV count; archs4
# / gtex (~18k / 21k genes) took ~20h / ~15h, recount2 (~6k genes) ~2.7h.
set -euo pipefail

DS="${1:?usage: 01_run_final_gls.sbatch <archs4|gtex|recount2> [--dry-run]}"
DRY="${2:-}"

CONDA_SH="${PHENOPLIER_CONDA_SH:-$HOME/miniforge3/etc/profile.d/conda.sh}"
ENV_NAME="${PHENOPLIER_ENV:-phenoplier-cli-neo}"
WS="${PHENOPLIER_HOME:-/pividori_lab/phenoplier_workspace}"
NAMESPACE="${PHENOPLIER_NAMESPACE:-clamp}"
TRAIT_FILTER="${CLAMP_TRAIT_FILTER:-biomedical}"
LV_PERCENTILE="${CLAMP_LV_PERCENTILE:-0.01}"

# shellcheck disable=SC1091
source "$CONDA_SH"
conda activate "$ENV_NAME"
export PATH="$CONDA_PREFIX/bin:$PATH"
export PHENOPLIER_HOME="$WS"
export PHENOPLIER_ROOT_DIR="$WS"

NJOBS="${SLURM_CPUS_PER_TASK:-16}"
echo "[start] $(date -Is)  host=$(hostname)  ds=$DS  njobs=$NJOBS  filter=$TRAIT_FILTER"
echo "[model] ${NAMESPACE}/final_${DS}_CLAMPfull_bp   workspace=$WS"
echo "[version] $(phenoplier -v 2>&1 | tail -1)"

time phenoplier shortcut gls \
  --model "${NAMESPACE}/final_${DS}_CLAMPfull_bp" \
  --lv-percentile "$LV_PERCENTILE" \
  --trait-filter "$TRAIT_FILTER" \
  --executor local \
  --n-jobs "$NJOBS" \
  ${DRY}
echo "[done] $(date -Is)"
