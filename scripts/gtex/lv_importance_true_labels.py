#!/usr/bin/env python3
# GTEx tissue feature importance analysis (true GTEx SMTSD labels + SHAP).
#
# Load GTEx LV data and tissue metadata
# Use true GTEx SMTSD tissue labels directly as one-vs-rest RF targets
# For each tissue, train a binary RF classifier (tissue vs. rest)
# Compute SHAP values to explain which LVs drive each tissue

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rpy2.robjects as ro

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import extract_B_matrix, extract_summary_matrix, run_tissue_rf_shap  # noqa: E402

readRDS = ro.r["readRDS"]

PARAM_DISTRIBUTIONS = {
    "n_estimators": [100, 300, 500, 1000],
    "max_depth": [None, 10, 20, 30],
    "min_samples_split": [2, 5, 10],
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clamp-full-rds", required=True)
    parser.add_argument("--gtex-meta", required=True)
    parser.add_argument("--min-samples", type=int, default=50)
    parser.add_argument("--global-seed", type=int, default=42)
    parser.add_argument("--min-rf-accuracy", type=float, default=0.9)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    np.random.seed(args.global_seed)

    MIN_SAMPLES = args.min_samples
    MIN_RF_ACCURACY = args.min_rf_accuracy

    OUTPUT_DIR = Path(args.out_dir)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"[gtex] Output directory: {OUTPUT_DIR}", flush=True)

    gtex_CLAMPfull_rds = readRDS(str(args.clamp_full_rds))
    gtex_CLAMPfull = extract_B_matrix(gtex_CLAMPfull_rds)
    lv_data = gtex_CLAMPfull

    gtex_CLAMPfull_summary = extract_summary_matrix(gtex_CLAMPfull_rds)
    gtex_CLAMPfull_summary["pathway"] = gtex_CLAMPfull_summary["pathway"].str.replace("C2CP_", "", regex=False)
    print(gtex_CLAMPfull_summary.shape, flush=True)

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
    print(gtex_meta.shape, flush=True)

    lv_matrix = lv_data.T

    meta_filtered = gtex_meta[gtex_meta["SAMPID"].isin(lv_matrix.index)].copy()
    meta_filtered = meta_filtered.set_index("SAMPID")
    meta_filtered = meta_filtered.loc[lv_matrix.index]

    tissue_counts = meta_filtered["SMTSD"].value_counts()
    valid_tissues = tissue_counts[tissue_counts >= MIN_SAMPLES].index
    mask = meta_filtered["SMTSD"].isin(valid_tissues)
    meta_filtered = meta_filtered[mask]
    lv_matrix = lv_matrix.loc[meta_filtered.index]

    tissue_labels = meta_filtered["SMTSD"]

    print(f"[gtex] Number of samples: {len(lv_matrix)}", flush=True)
    print(f"[gtex] Number of LVs: {lv_matrix.shape[1]}", flush=True)
    print(f"[gtex] Number of unique tissues (>= {MIN_SAMPLES} samples): {tissue_labels.nunique()}", flush=True)

    true_tissue_labels = tissue_labels.copy()
    print("[gtex] GTEx true tissue label distribution:", flush=True)
    print(true_tissue_labels.value_counts(), flush=True)

    tissues = sorted(true_tissue_labels.unique())
    print(f"[gtex] GTEx true-label tissues: {tissues}", flush=True)

    def extra_accuracy_fields(tissue: str) -> dict:
        n_pos = int((true_tissue_labels == tissue).sum())
        n_neg = int((true_tissue_labels != tissue).sum())
        return {
            "N_Positive": n_pos,
            "N_Negative": n_neg,
            "Label_Source": "GTEx_SMTSD",
        }

    run_tissue_rf_shap(
        lv_matrix,
        tissues,
        label_fn=lambda tissue: (true_tissue_labels == tissue),
        output_dir=OUTPUT_DIR,
        global_seed=args.global_seed,
        min_rf_accuracy=MIN_RF_ACCURACY,
        param_distributions=PARAM_DISTRIBUTIONS,
        extra_accuracy_fields_fn=extra_accuracy_fields,
    )


if __name__ == "__main__":
    main()
