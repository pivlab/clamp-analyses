#!/usr/bin/env python3
"""GTEx subtissue-recovery donor-grouped fold assignment (runs once).

Fixes the valid sample universe (SMTS categories with >= min-samples samples,
same filter as the main true-labels pipeline) and assigns every sample to one
of n-folds donor-grouped folds, stratified on the 27-class SMTS label so no
downstream one-vs-rest RF's negative class is skewed by an uneven fold split.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rpy2.robjects as ro
from sklearn.model_selection import StratifiedGroupKFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import derive_donor_id, extract_B_matrix  # noqa: E402

readRDS = ro.r["readRDS"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clamp-full-rds", required=True)
    parser.add_argument("--gtex-meta", required=True)
    parser.add_argument("--min-samples", type=int, default=50)
    parser.add_argument("--n-folds", type=int, default=5)
    parser.add_argument("--global-seed", type=int, default=42)
    parser.add_argument("--out-membership", required=True)
    parser.add_argument("--out-summary", required=True)
    args = parser.parse_args()

    gtex_CLAMPfull_rds = readRDS(str(args.clamp_full_rds))
    lv_matrix = extract_B_matrix(gtex_CLAMPfull_rds).T

    gtex_meta = pd.read_csv(
        args.gtex_meta,
        sep="\t",
        header=0,
        dtype=str,
        quoting=csv.QUOTE_NONE,
        engine="python",
        comment=None,
        keep_default_na=False,
        on_bad_lines="warn",
    )
    print(f"[gtex] {gtex_meta.shape}", flush=True)

    meta_filtered = gtex_meta[gtex_meta["SAMPID"].isin(lv_matrix.index)].copy()
    meta_filtered = meta_filtered.set_index("SAMPID")
    meta_filtered = meta_filtered.loc[lv_matrix.index]
    meta_filtered.index.name = "SAMPID"  # .loc against lv_matrix.index (unnamed) clears the index name

    # Broad-tissue universe filter: always SMTS, regardless of what the main
    # true-labels pipeline uses, since this analysis is about subtissue
    # recoverability from a tissue-blind (SMTS) RF.
    tissue_counts = meta_filtered["SMTS"].value_counts()
    valid_tissues = tissue_counts[tissue_counts >= args.min_samples].index
    meta_filtered = meta_filtered[meta_filtered["SMTS"].isin(valid_tissues)]
    meta_filtered = meta_filtered.reset_index()
    meta_filtered["donor_id"] = derive_donor_id(meta_filtered["SAMPID"])

    print(f"[gtex] Valid samples: {len(meta_filtered)}", flush=True)
    print(f"[gtex] Valid SMTS tissues (>= {args.min_samples} samples): {len(valid_tissues)}", flush=True)
    print(f"[gtex] Unique donors: {meta_filtered['donor_id'].nunique()}", flush=True)

    splitter = StratifiedGroupKFold(n_splits=args.n_folds, shuffle=True, random_state=args.global_seed)
    fold_of = np.zeros(len(meta_filtered), dtype=int)
    for fold_idx, (_, test_idx) in enumerate(
        splitter.split(meta_filtered, meta_filtered["SMTS"], groups=meta_filtered["donor_id"])
    ):
        fold_of[test_idx] = fold_idx + 1

    meta_filtered["fold"] = fold_of
    meta_filtered["split_seed"] = args.global_seed

    membership = meta_filtered[["SAMPID", "donor_id", "SMTS", "SMTSD", "fold", "split_seed"]].copy()

    donor_fold_counts = membership.groupby("donor_id")["fold"].nunique()
    assert (donor_fold_counts == 1).all(), "a donor spans multiple folds"
    assert set(membership["fold"].unique()) == set(range(1, args.n_folds + 1)), "at least one fold is empty"

    out_membership = Path(args.out_membership)
    out_summary = Path(args.out_summary)
    out_membership.parent.mkdir(parents=True, exist_ok=True)
    out_summary.parent.mkdir(parents=True, exist_ok=True)

    membership.to_csv(out_membership, sep="\t", index=False)

    per_tissue_summary = (
        membership.groupby(["fold", "SMTS"])
        .agg(n_samples=("SAMPID", "size"), n_donors=("donor_id", "nunique"))
        .reset_index()
    )
    totals = (
        membership.groupby("fold")
        .agg(n_samples=("SAMPID", "size"), n_donors=("donor_id", "nunique"))
        .reset_index()
    )
    totals["SMTS"] = "TOTAL"
    fold_summary = pd.concat([per_tissue_summary, totals], ignore_index=True)
    fold_summary = fold_summary.sort_values(["fold", "SMTS"]).reset_index(drop=True)
    fold_summary.to_csv(out_summary, sep="\t", index=False)

    print(f"[gtex] Wrote fold membership to {out_membership}", flush=True)
    print(f"[gtex] Wrote fold summary to {out_summary}", flush=True)


if __name__ == "__main__":
    main()
