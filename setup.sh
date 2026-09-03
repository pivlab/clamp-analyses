#!/usr/bin/env bash
# Bootstraps the local dev machine for this repo: installs Conda if missing,
# then creates (idempotently) the three named Conda environments every
# Snakemake rule's `conda:` directive assumes already exist:
#   - clamp-analyses  (envs/clamp-analyses.lock)  -- default CPU env
#   - gpu-kmeans      (envs/gpu-kmeans.lock)       -- GPU clustering env
#   - snakemake       (envs/snakemake.yaml)        -- pipeline orchestrator
# Re-running this script is safe: existing envs are left untouched and the
# external CLAMP package is verified against the repository's pinned revision.
#
# Lock files are platform-pinned (linux-64) -- this script currently only
# supports Linux x86_64.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_ROOT"

MINICONDA_PREFIX="${MINICONDA_PREFIX:-$HOME/miniconda3}"

log() { printf '\n[setup.sh] %s\n' "$1"; }

# ------------------------------------------------------------------
# 1. Ensure conda is available (install Miniconda if not)
# ------------------------------------------------------------------
if command -v conda >/dev/null 2>&1; then
    log "conda already installed: $(command -v conda)"
    CONDA_BASE="$(conda info --base)"
else
    log "conda not found, installing Miniconda to $MINICONDA_PREFIX"
    case "$(uname -s)-$(uname -m)" in
        Linux-x86_64) installer="Miniconda3-latest-Linux-x86_64.sh" ;;
        *)
            echo "[setup.sh] unsupported platform $(uname -s)-$(uname -m)." >&2
            echo "[setup.sh] this repo's env lock files are pinned to linux-64;" >&2
            echo "[setup.sh] install Conda/Mambaforge yourself and re-run this script." >&2
            exit 1
            ;;
    esac
    tmp_installer="$(mktemp -d)/miniconda.sh"
    curl -fsSL "https://repo.anaconda.com/miniconda/${installer}" -o "$tmp_installer"
    bash "$tmp_installer" -b -p "$MINICONDA_PREFIX"
    rm -f "$tmp_installer"
    CONDA_BASE="$MINICONDA_PREFIX"
fi

# shellcheck disable=SC1091
source "$CONDA_BASE/etc/profile.d/conda.sh"

# ------------------------------------------------------------------
# 2. Create the three named envs, skipping any that already exist
# ------------------------------------------------------------------
env_exists() {
    conda env list | awk '{print $1}' | grep -qx "$1"
}

create_env_from_lock() {
    local name="$1" lock="$2"
    if env_exists "$name"; then
        log "env '$name' already exists, skipping (remove it manually with 'conda env remove -n $name' to recreate)"
    else
        log "creating env '$name' from $lock"
        conda create --name "$name" --file "$lock" -y
    fi
}

create_env_from_yaml() {
    local name="$1" yaml="$2"
    if env_exists "$name"; then
        log "env '$name' already exists, skipping (remove it manually with 'conda env remove -n $name' to recreate)"
    else
        log "creating env '$name' from $yaml"
        conda env create -n "$name" -f "$yaml"
    fi
}

create_env_from_lock clamp-analyses envs/clamp-analyses.lock
create_env_from_lock gpu-kmeans envs/gpu-kmeans.lock
create_env_from_yaml snakemake envs/snakemake.yaml

# ------------------------------------------------------------------
# 3. Install/verify the pinned external CLAMP R package
# ------------------------------------------------------------------
log "WARNING: using the repository-pinned CLAMP revision; do not update CLAMP independently"
conda run -n clamp-analyses Rscript scripts/install_clamp.R

# ------------------------------------------------------------------
# 4. Done
# ------------------------------------------------------------------
log "setup complete. Next steps:"
cat <<'EOF'

  conda activate snakemake
  snakemake --cores 4 --use-conda --snakefile workflow/Snakefile --touch all

To verify everything (envs, Snakemake, data files), run:
  nbs/00_setup/00_check_setup.ipynb
EOF
