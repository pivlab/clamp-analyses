#!/usr/bin/env bash
# Concatenate every per-model GLS summary, then drop excluded UK Biobank
# phenotype categories. Runs in the phenoplier-cli-neo conda env: the
# category filter (filter_phenotypes.py) needs phenoplier's own
# phenoplier.entity.Trait.category, not this repo's clamp-analyses env.
set -euo pipefail

root="$1"
raw_out="$2"
final_out="$3"
excluded_log="$4"
conda_env="$5"
shift 5
exclude_terms=("$@")

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$conda_env"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python "$script_dir/aggregate_traits.py" --root "$root" --out "$raw_out"
python "$script_dir/filter_phenotypes.py" \
  --input "$raw_out" --output "$final_out" --excluded-log "$excluded_log" \
  --exclude "${exclude_terms[@]}"
