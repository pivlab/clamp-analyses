#!/usr/bin/env bash
#
# Move the pre-existing ARCHS4 results into the layout the Snakemake rules
# declare, so `snakemake --touch archs4_precomputed` can adopt them instead of
# recomputing 2.2 TB.
#
# Everything here is a rename within one filesystem, so each move is O(1) and
# the volume never needs room for a second copy.  Nothing is ever deleted:
# tracks we no longer analyse are parked under output/_deprecated/ so they stay
# recoverable.
#
# Usage:
#   scripts/archs4/migrate_outputs.sh            # dry run, prints the plan
#   scripts/archs4/migrate_outputs.sh --apply    # perform the moves
#
# Re-running after a successful --apply is a no-op: a pair whose source is gone
# and whose destination exists is reported as already migrated.

set -euo pipefail

APPLY=0
for arg in "$@"; do
    case "$arg" in
        --apply) APPLY=1 ;;
        -h|--help) sed -n '2,20p' "$0"; exit 0 ;;
        *) echo "unknown argument: $arg" >&2; exit 2 ;;
    esac
done

REPO_ROOT="$(git -C "$(dirname "${BASH_SOURCE[0]}")" rev-parse --show-toplevel)"
cd "$REPO_ROOT"

MB="output/01_model_building"
BIO="output/03_model_biology"
PROD="$MB/02_archs4"
NEWBIO="$BIO/02_archs4"
DEPRECATED="output/_deprecated"
LOG="output/_migration_log.txt"

# ---------------------------------------------------------------------------
# The migration table.  One "source|destination" pair per line; order matters
# only in that the final model must land before anything referring to it.
# ---------------------------------------------------------------------------
MOVES=(
    # -- the shared preprocessing output: the filtered, z-scored FBM that every
    #    ARCHS4 model is fit on, plus the gene/sample metadata describing it.
    "$MB/04_archs4/01_archs4_preprocess|$PROD/00_preprocess"

    # -- the canonical full model (K=1728).  Two container levels collapse; the
    #    hall_coverage_rs100_seed_{1,2,3} leaves inside keep their names.  Each
    #    leaf is self-contained (its own svd.rds, CLAMP_K.rds, CLAMPbase.rds and
    #    CLAMPfull_hall.rds), which is why 02_archs4_svd / 03_archs4_CLAMPbase /
    #    04_archs4_CLAMPfull below are superseded rather than upstream of it.
    "$MB/04_archs4/06_bp_coverage_rshall/06_bp_coverage_hall_rs_100|$PROD/01_final_model"

    # -- subsampling sweeps we keep (by study) ------------------------------
    "$MB/04_archs4/07_bp_coverage_study|$PROD/02_coverage_study"
    "$MB/04_archs4/08_saturation_study|$PROD/03_saturation_study"

    # -- reconstruction error and the LINCS projection ----------------------
    "$MB/06_reconstruction_error|$PROD/04_reconstruction_error"
    "$MB/07_drug_perturbation/lincs|$PROD/05_drug_disease/lincs"

    # -- biology: coverage (ORA) --------------------------------------------
    "$BIO/00_archs4/06_coverage_study/gpu_ora_bp|$NEWBIO/00_coverage/gpu_ora_bp"
    "$BIO/00_archs4/06_coverage_study/gpu_ora_msigdb|$NEWBIO/00_coverage/gpu_ora_msigdb"

    # -- biology: cross-dataset ARCHS4 vs GTEx vs recount2.  These were filed
    #    under the saturation tree, which is where they do not belong.
    "$BIO/00_archs4/08_saturation_study/gpu_ora_bp/GTEx|$NEWBIO/00_coverage/cross_dataset/gpu_ora_bp/GTEx"
    "$BIO/00_archs4/08_saturation_study/gpu_ora_bp/Recount2|$NEWBIO/00_coverage/cross_dataset/gpu_ora_bp/Recount2"
    "$BIO/00_archs4/08_saturation_study/gpu_ora_bp/results_external_bp_gpu_ora.csv|$NEWBIO/00_coverage/cross_dataset/gpu_ora_bp/results_external_bp_gpu_ora.csv"
    "$BIO/00_archs4/08_saturation_study/gpu_ora_shared_universe|$NEWBIO/00_coverage/cross_dataset/shared_universe"

    # -- biology: saturation (ORA).  Must follow the cross-dataset moves
    #    above, which carve GTEx/Recount2 out of these same directories.
    "$BIO/00_archs4/08_saturation_study/gpu_ora_bp|$NEWBIO/01_saturation/gpu_ora_bp"
    "$BIO/00_archs4/08_saturation_study/gpu_ora_msigdb|$NEWBIO/01_saturation/gpu_ora_msigdb"

    # -- biology: the remaining analyses we keep ----------------------------
    "$BIO/00_archs4/02_drug_disease_associations|$NEWBIO/02_drug_diseases_associations"
    "$BIO/00_archs4/03_tissue_projections|$NEWBIO/03_projections"
    "$BIO/00_archs4/01_CRISPRCas9|$NEWBIO/04_crispercas"
)

# Tracks that stay on disk but leave the analysis.  Random-sample coverage and
# saturation, the ORA implementations we did not keep (decoupler, blitzgsea),
# the robustness study, and the plot directories full of committed PNG/PDF --
# the new notebooks must not emit images.
DEPRECATE=(
    "$MB/04_archs4/00_archs4_download"
    # The superseded base chain.  Its CLAMPfull step wrote CLAMPbase's loadings
    # into the full result, so 04_archs4_CLAMPfull/CLAMPfull/Z.csv is byte-for-
    # byte CLAMPbase/Z.csv.  Nothing downstream ever read it -- the shipped model
    # is 01_final_model, fit by the coverage sweep, whose code is correct.
    "$MB/04_archs4/02_archs4_svd"
    "$MB/04_archs4/03_archs4_CLAMPbase"
    "$MB/04_archs4/04_archs4_CLAMPfull"
    "$MB/04_archs4/06_bp_coverage_rshall"
    "$MB/04_archs4/07_saturation"
    "$BIO/00_archs4/00_pathway_coverage"
    "$BIO/00_archs4/02_pathway_coverage_hall_study"
    "$BIO/00_archs4/02_tissue_archs4"
    "$BIO/00_archs4/03_biology_validation_RF"
    "$BIO/00_archs4/03_pathway_coverage_hall_saturation_study"
    "$BIO/00_archs4/04_pseudobulk"
    "$BIO/00_archs4/05_coverage_random"
    "$BIO/00_archs4/07_saturation_random"
    "$BIO/00_archs4/10_archs4_robutsness"
    "$BIO/00_archs4/_blitzgsea_100pct_anchor"
    "$BIO/00_archs4/06_coverage_study"
    "$BIO/00_archs4/08_saturation_study"
)

# Files inside the reconstruction-error directory that belong to the
# random-sample track or are rendered figures.  These are resolved after that
# directory has moved; on a dry run it has not moved yet, so the patterns are
# matched against whichever of the two locations currently exists.
RECON_DIR_SRC="$MB/06_reconstruction_error"
RECON_DIR_DST="$PROD/04_reconstruction_error"
DEPRECATE_RECON_PATTERNS=( "*random*" "*.pdf" )

# ---------------------------------------------------------------------------

n_move=0 n_skip=0 n_missing=0 n_error=0

log() {
    if (( APPLY )); then printf '%s\n' "$*" >>"$LOG"; fi
}

# Device id of a path, or of its nearest existing ancestor when it does not
# exist yet -- a destination is normally created by this script.
device_of() {
    local p="$1"
    while [[ ! -e "$p" && "$p" != "." && "$p" != "/" ]]; do p="$(dirname "$p")"; done
    stat -c %d "$p"
}

plan_move() {
    local src="$1" dst="$2"

    if [[ ! -e "$src" ]]; then
        if [[ -e "$dst" ]]; then
            echo "  already migrated: $dst"
            (( ++n_skip ))
        else
            echo "  MISSING (neither side exists): $src"
            (( ++n_missing ))
        fi
        return
    fi

    if [[ -e "$dst" ]]; then
        echo "  ERROR destination exists, refusing to merge: $dst"
        (( ++n_error ))
        return
    fi

    local dsrc ddst
    dsrc="$(device_of "$src")"
    ddst="$(device_of "$dst")"
    if [[ "$dsrc" != "$ddst" ]]; then
        echo "  ERROR cross-device move ($dsrc -> $ddst); this would copy, not rename:"
        echo "        $src -> $dst"
        (( ++n_error ))
        return
    fi

    echo "  mv $src -> $dst"
    (( ++n_move ))
    if (( APPLY )); then
        mkdir -p "$(dirname "$dst")"
        mv "$src" "$dst"
        log "$(date -Is) mv $src -> $dst"
    fi
}

echo "ARCHS4 output migration -- $( ((APPLY)) && echo APPLY || echo 'DRY RUN (pass --apply to execute)')"
echo
echo "Analyses we keep:"
for pair in "${MOVES[@]}"; do
    plan_move "${pair%%|*}" "${pair##*|}"
done

echo
echo "Tracks moved aside to $DEPRECATED (kept, not deleted):"
for src in "${DEPRECATE[@]}"; do
    plan_move "$src" "$DEPRECATED/${src#output/}"
done

echo
echo "Random-track and figure files inside the reconstruction-error directory:"
# On --apply the directory has already moved; on a dry run it has not.  Preview
# against whichever location holds the files right now, but always deprecate to
# one fixed path so the outcome does not depend on when this ran.
recon_dir="$RECON_DIR_DST"
[[ -d "$recon_dir" ]] || recon_dir="$RECON_DIR_SRC"
shopt -s nullglob
for pattern in "${DEPRECATE_RECON_PATTERNS[@]}"; do
    for f in "$recon_dir"/$pattern; do
        plan_move "$f" "$DEPRECATED/01_model_building/06_reconstruction_error/$(basename "$f")"
    done
done
shopt -u nullglob

echo
echo "moves: $n_move   already migrated: $n_skip   missing: $n_missing   errors: $n_error"

if (( n_error )); then
    echo "Refusing to proceed while errors remain." >&2
    exit 1
fi

if (( APPLY )); then
    # The now-empty shells of the old layout.  rmdir, never rm -r: it fails
    # loudly if anything unexpected is still in there.
    for d in "$MB/04_archs4" "$BIO/00_archs4" "$MB/07_drug_perturbation"; do
        [[ -d "$d" ]] && rmdir "$d" 2>/dev/null && echo "removed empty $d"
    done
    echo "Done. Log: $LOG"
else
    echo "Nothing was changed. Re-run with --apply to perform these moves."
fi
