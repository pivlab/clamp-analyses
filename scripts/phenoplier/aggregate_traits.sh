#!/usr/bin/env bash
# Concatenate per-model GLS summaries into one long table, inside the
# phenoplier-cli conda env for a consistent pandas.
#
#   aggregate_traits.sh <conda-env> --out <long.csv> \
#       --clampfull-bp <summaries...> --clampbase <summaries...>
#
# Non-biomedical traits (job codes, medication codes, diet items, ...) are
# excluded by phenoplier-cli itself at run time -- see phenoplier.yaml's
# trait_filter and run_gls.sh's --trait-filter -- so they never reach these
# summaries and there is nothing to filter here. Each model's exclusion log
# travels next to its summary as <name>_trait_filter_excluded.tsv.
set -euo pipefail

conda_env="$1"
shift

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$conda_env"

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python "$script_dir/aggregate_traits.py" "$@"
