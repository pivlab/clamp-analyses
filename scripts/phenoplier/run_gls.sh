#!/usr/bin/env bash
# Run `phenoplier shortcut gls` for one CLAMP model and copy its result back
# into this repo's output tree.
#
# phenoplier-cli lives in its own conda env (see setup_env.sh) and its
# `shortcut gls` names the project directory itself (date-prefixed stem, see
# its docs/tutorials/tutorial-custom-model.md); rather than reproducing that
# naming logic here, this script passes an explicit --name and then finds the
# resulting project directory by that name, so a change to phenoplier-cli's
# naming scheme doesn't silently break this wrapper.
set -euo pipefail

rds="$1"
model_key="$2"
conda_env="$3"
namespace="$4"
lv_percentile="$5"
executor="$6"
cluster="$7"
n_jobs="$8"
summary_out="$9"

name="$(printf '%s' "$model_key" | tr '/' '_')"
workspace="${PHENOPLIER_HOME:-$HOME/phenoplier}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$conda_env"

cluster_args=()
if [[ -n "$cluster" ]]; then
  cluster_args=(--cluster "$cluster")
fi

phenoplier shortcut gls \
  --input "$rds" \
  --name "$name" \
  --model-namespace "$namespace" \
  --lv-percentile "$lv_percentile" \
  --executor "$executor" \
  --n-jobs "$n_jobs" \
  "${cluster_args[@]}"

project_dir="$(
  find "$workspace/projects" -maxdepth 1 -type d -name "*${name}*" \
    -printf '%T@ %p\n' | sort -rn | head -n1 | cut -d' ' -f2-
)"
if [[ -z "$project_dir" ]]; then
  echo "No project directory matching '*${name}*' found under $workspace/projects" >&2
  exit 1
fi

mkdir -p "$(dirname "$summary_out")"
cp "$project_dir/results/gls/gls-summary-phenomexcan.tsv.gz" "$summary_out"
cp "$project_dir/results/gls/gls-summary-phenomexcan.pkl.gz" "$(dirname "$summary_out")/" \
  2>/dev/null || true
