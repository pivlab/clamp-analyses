#!/usr/bin/env bash
# Build a one-file HDF5 composite study store for one CLAMP model: the model's
# .rds (Z matrix, hyperparameters, provenance) together with that model's GLS
# results for one cohort. Wraps `phenoplier store build` (phenoplier-cli #70).
# Driven by workflow/rules/phenoplier.smk; the project is the one run_gls.sh
# created under the same --name, so the two stay in step.
#
# Built atomically: to a .tmp, verified to open with h5py, then renamed, so a
# failed build never leaves a partial .h5 for Snakemake to mistake for done.
set -euo pipefail

usage() {
  cat >&2 <<USAGE
usage: $0 --rds <model.rds> --name <project> --store-out <out.h5> [--conda-env NAME] [--require-version X.Y.Z] [--cohort NAME]
Workspace: \$PHENOPLIER_HOME (default \$HOME/phenoplier).
USAGE
  exit 2
}

conda_env=phenoplier-cli-neo
require_version=""
cohort=phenomexcan_rapid_gwas
rds=""
name=""
store_out=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --rds) rds="$2"; shift 2 ;;
    --name) name="$2"; shift 2 ;;
    --store-out) store_out="$2"; shift 2 ;;
    --conda-env) conda_env="$2"; shift 2 ;;
    --require-version) require_version="$2"; shift 2 ;;
    --cohort) cohort="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "unknown argument: $1" >&2; usage ;;
  esac
done
[[ -n "$rds" && -n "$name" && -n "$store_out" ]] || usage
[[ -f "$rds" ]] || { echo "model not found: $rds" >&2; exit 1; }

workspace="${PHENOPLIER_HOME:-$HOME/phenoplier}"
export PHENOPLIER_HOME="$workspace" PHENOPLIER_ROOT_DIR="$workspace"

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$conda_env"
# Reading the .rds via `store build --clamp-rds` needs R on PATH too (rpy2).
export PATH="${CONDA_PREFIX:-}/bin:${PATH}"

installed="$(phenoplier --version 2>&1 | tail -1 | sed -E 's/.*v([0-9][0-9A-Za-z.+-]*).*/\1/')"
if [[ -n "$require_version" && "$installed" != "$require_version" ]]; then
  echo "[ERROR] phenoplier-cli $installed in env '$conda_env' but $require_version is required" \
    "(workflow/config/phenoplier.yaml: version)" >&2
  exit 1
fi

project="$workspace/projects/$name"
# Per-phenotype GLS results live under the cohort-named subdir; the combined
# gls-summary lives one level up (store build auto-discovers it there).
gls_dir="$project/results/gls/phenoplier/${cohort}"
if [[ ! -d "$gls_dir" ]]; then
  echo "[ERROR] $name: no GLS results dir at $gls_dir" >&2
  exit 1
fi

mkdir -p "$(dirname "$store_out")"
tmp="${store_out}.tmp.$$"
rm -f "$tmp"
echo "[store] $(date -Is)  $name -> $store_out"
phenoplier store build \
  --output "$tmp" \
  --cohort "$cohort" \
  --clamp-rds "$rds" \
  --gls-dir "$gls_dir" \
  --method-name CLAMP
python -c "import h5py, sys; h5py.File(sys.argv[1], 'r').close()" "$tmp"
mv -f "$tmp" "$store_out"
echo "[done] $(date -Is)  $name"
