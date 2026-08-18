#!/usr/bin/env bash
# One-time setup of the isolated conda env phenoplier-cli runs in on this
# machine/cluster. Not called by Snakemake -- run once per machine, same as
# any other manual environment setup in this repo (envs/*.yaml).
#
# Delegates to phenoplier-cli's own installer rather than reimplementing it;
# see https://github.com/pivlab/phenoplier-cli#installation.
set -euo pipefail

# Pinned to the open fix for the GLS standard-error collapse (PR #85) rather
# than `main`: that bug produces artifactual ~0 p-values specifically in the
# hall_coverage_rs*-style sub-sampling models this integration runs against
# (see workflow/rules/phenoplier.smk). Switch back to `main` once #85 merges.
REF="${PHENOPLIER_REF:-fix/gls-se-collapse-guard-reopen}"
CLONE_DIR="${PHENOPLIER_CLONE_DIR:-$HOME/phenoplier-cli}"
TARGET_ENV="${PHENOPLIER_TARGET_ENV:-phenoplier-cli-neo}"

# install_phenoplier_cli.sh clones/updates $CLONE_DIR itself if needed --
# see https://github.com/pivlab/phenoplier-cli/blob/main/scripts/install_phenoplier_cli.sh
curl -fsSL "https://raw.githubusercontent.com/pivlab/phenoplier-cli/main/scripts/install_phenoplier_cli.sh" \
  | bash -s -- --dir "$CLONE_DIR" --env "$TARGET_ENV" --ref "$REF" --extras slurm "$@"

cat <<EOF

Next:
  1. Sanity-check the install itself (no full data bundle needed):
       bash scripts/phenoplier/verify_install.sh $TARGET_ENV

  2. Wire the reference data bundle into the workspace:
       conda activate $TARGET_ENV
       phenoplier workspace init
       phenoplier workspace link /path/to/phenoplier_full_data
         (Alpine:    /pl/active/pivlab/projects/hzhang/data/phenoplier_full_data.tar.gz)
         (server_cu: /pividori_lab/data/phenoplier_full_data.tar.gz)

  3. Smoke-test this repo's integration on a small (K=86 LV) model:
       snakemake --profile workflow/profiles/local \\
         output/03_model_biology/00_phenoplier/archs4_saturation/rs1/k86/seed1/CLAMPbase/gls-summary-phenomexcan.tsv.gz
EOF
