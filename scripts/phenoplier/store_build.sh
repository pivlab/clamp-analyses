#!/usr/bin/env bash
# Build a one-file HDF5 composite study store for one CLAMP model: the model's
# .rds (Z matrix, hyperparameters, provenance) together with that model's GLS
# results for the phenomexcan cohort. Wraps `phenoplier store build`
# (phenoplier-cli #70).
#
# Like run_gls.sh, this finds the phenoplier-cli project directory by the name
# derived from the model_key rather than reproducing the CLI's naming logic, so
# the two stay in step.
set -euo pipefail

rds="$1"
model_key="$2"
conda_env="$3"
cohort="$4"
store_out="$5"

name="$(printf '%s' "$model_key" | tr '/' '_')"
workspace="${PHENOPLIER_HOME:-$HOME/phenoplier}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$conda_env"
# `phenoplier` shells out to bare executables; keep the env bin on PATH. Reading
# the .rds via `store build --clamp-rds` needs R on PATH too (rpy2 -> R RHOME).
export PATH="${CONDA_PREFIX:-}/bin:${PATH}"

project_dir="$(
  find "$workspace/projects" -maxdepth 1 -type d -name "*${name}*" \
    -printf '%T@ %p\n' | sort -rn | head -n1 | cut -d' ' -f2-
)"
if [[ -z "$project_dir" ]]; then
  echo "No project directory matching '*${name}*' found under $workspace/projects" >&2
  exit 1
fi

# Per-phenotype GLS results live under the cohort-named subdir; the combined
# gls-summary lives one level up (store build auto-discovers it there).
gls_dir="$project_dir/results/gls/phenoplier/${cohort}"
if [[ ! -d "$gls_dir" ]]; then
  echo "No GLS results dir at $gls_dir" >&2
  exit 1
fi

mkdir -p "$(dirname "$store_out")"

phenoplier store build \
  --output "$store_out" \
  --cohort "$cohort" \
  --clamp-rds "$rds" \
  --gls-dir "$gls_dir" \
  --method-name CLAMP
