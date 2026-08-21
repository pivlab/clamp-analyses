#!/usr/bin/env bash
# Concatenate every per-model GLS summary into one long table.
#
# Non-biomedical traits (job codes, medication codes, diet items, ...) are
# excluded by phenoplier-cli itself at run time -- see phenoplier.yaml's
# trait_filter and run_gls.sh's --trait-filter -- so they never reach these
# summaries and there is nothing to filter here. Each model's exclusion log
# travels next to its summary as trait_filter_excluded.tsv.
#
# Still runs in the phenoplier-cli-neo conda env for a consistent pandas.
set -euo pipefail

root="$1"
final_out="$2"
conda_env="$3"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$conda_env"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python "$script_dir/aggregate_traits.py" --root "$root" --out "$final_out"
