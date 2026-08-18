#!/usr/bin/env bash
# Verify a phenoplier-cli install is healthy, independent of this repo's own
# models: runs phenoplier-cli's own test suite against its small `test_data`
# bundle, mirroring its CI (.github/workflows/main.yml) rather than the full
# ~100+ GB reference data bundle. Run once after setup_env.sh, before trying
# a real GLS run against one of this repo's CLAMP models.
set -euo pipefail

conda_env="${1:-phenoplier-cli-neo}"
clone_dir="${PHENOPLIER_CLONE_DIR:-$HOME/phenoplier-cli}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$conda_env"

python -m pip show pytest >/dev/null 2>&1 || python -m pip install pytest pytest-cov pytest-env

phenoplier workspace init
phenoplier workspace download test_data --parallel

(cd "$clone_dir" && PYTHONPATH=. pytest -rs -m "not requires_full_data" test/)
