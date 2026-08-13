#!/usr/bin/env bash
#
# Move code up to the lab server and ARCHS4 results back down.
#
#   scripts/archs4/sync_server.sh push    # code, configs, notebooks -> server
#   scripts/archs4/sync_server.sh pull    # ARCHS4 outputs -> workstation
#
# The ARCHS4 CLAMP fits need ~500 GB of RAM, which no workstation here has, so
# they run on server_cu under workflow/profiles/slurm. The GPU ORA sweeps and
# every plotting notebook run locally. This script is the seam between the two.
#
# Add --dry-run to preview either direction (it is forwarded to rsync).

set -euo pipefail

REMOTE_HOST="${ARCHS4_REMOTE_HOST:-server_cu}"
REMOTE_ROOT="${ARCHS4_REMOTE_ROOT:-/pividori_lab/marc_projects/clamp-analyses}"

REPO_ROOT="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "$REPO_ROOT"

# Results live only where they were computed: never let a push overwrite the
# server's outputs, and never let a pull clobber local ones outside ARCHS4.
PUSH_EXCLUDES=(
    --exclude 'output/'
    --exclude 'data/'
    --exclude '.git/'
    --exclude '.snakemake/'
    --exclude '__pycache__/'
    --exclude '.ipynb_checkpoints/'
    --exclude '*.executed.ipynb'
)

# The ARCHS4 trees the server produces. Deliberately narrow: pulling all of
# output/ would drag down other datasets' results and the _deprecated tree.
PULL_PATHS=(
    "output/01_model_building/02_archs4/"
)

direction="${1:-}"
shift || true

case "$direction" in
    push)
        echo "push: $REPO_ROOT -> $REMOTE_HOST:$REMOTE_ROOT"
        rsync -avz --delete-after "${PUSH_EXCLUDES[@]}" "$@" \
            ./ "$REMOTE_HOST:$REMOTE_ROOT/"
        ;;
    pull)
        for path in "${PULL_PATHS[@]}"; do
            echo "pull: $REMOTE_HOST:$REMOTE_ROOT/$path -> $path"
            mkdir -p "$path"
            # No --delete: the server holds a subset of what we have locally,
            # and a stray delete here would cost days of compute to rebuild.
            rsync -avz "$@" "$REMOTE_HOST:$REMOTE_ROOT/$path" "$path"
        done
        echo
        echo "Now adopt the pulled results so Snakemake does not refit them:"
        echo "  snakemake --profile workflow/profiles/local --touch archs4_precomputed"
        ;;
    *)
        echo "usage: $(basename "$0") {push|pull} [rsync options]" >&2
        exit 2
        ;;
esac
