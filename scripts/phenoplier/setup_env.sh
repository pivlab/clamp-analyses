#!/usr/bin/env bash
# One-time setup of the isolated conda env phenoplier-cli runs in on this
# machine/cluster. Not called by Snakemake -- run once per machine, same as
# any other manual environment setup in this repo (envs/*.yaml).
#
# Delegates to phenoplier-cli's own installer rather than reimplementing it;
# see https://github.com/pivlab/phenoplier-cli#installation.
set -euo pipefail

# Pinned to a release tag so every machine installs the same code. v0.5.1 is
# the first release carrying both pieces this integration depends on:
#   * the GLS standard-error collapse fix (phenoplier-cli #85), which removes
#     the artifactual ~0 p-values that hit the sub-sampling models this runs
#     against (see workflow/rules/phenoplier.smk);
#   * the trait filter (#99 + the #100 settings fix), used via
#     `--trait-filter` in run_gls.sh.
REF="${PHENOPLIER_REF:-v0.5.1}"
CLONE_DIR="${PHENOPLIER_CLONE_DIR:-$HOME/phenoplier-cli}"
TARGET_ENV="${PHENOPLIER_TARGET_ENV:-phenoplier-cli-neo}"

# install_phenoplier_cli.sh clones/updates $CLONE_DIR itself if needed --
# see https://github.com/pivlab/phenoplier-cli/blob/main/scripts/install_phenoplier_cli.sh
curl -fsSL "https://raw.githubusercontent.com/pivlab/phenoplier-cli/main/scripts/install_phenoplier_cli.sh" \
  | bash -s -- --dir "$CLONE_DIR" --env "$TARGET_ENV" --ref "$REF" --extras slurm "$@"

cat <<EOF

Next: wire the reference data bundle into the workspace once:
  conda activate $TARGET_ENV
  phenoplier workspace init
  phenoplier workspace link /path/to/phenoplier_full_data
    (Alpine:    /pl/active/pivlab/projects/hzhang/data/phenoplier_full_data.tar.gz)
    (server_cu: /pividori_lab/data/phenoplier_full_data.tar.gz)
EOF
