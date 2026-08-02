#!/usr/bin/env python3
"""GTEx tissue feature importance analysis (K-means binary + SHAP).

1. Load GTEx LV data, tissue metadata, and the saved K-means model
2. Derive binary tissue labels from K-means cluster assignments (majority-vote
   tissue per cluster)
3. For each tissue, train a binary RF classifier (tissue vs. rest) using
   k-means-derived labels
4. Compute SHAP values to explain which LVs drive each tissue

Training half of nbs/03_model_biology/01_gtex/02_rf_kmeans/00_LV_importance_kmeans.ipynb;
the notebook itself only reads this script's output TSVs and plots the LV
correlation heatmap.
"""

from __future__ import annotations

import argparse
import csv
import pickle
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri
from sklearn.metrics import adjusted_rand_score

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
    parser.add_argument("--kmeans-model", required=True)
    parser.add_argument("--min-samples", type=int, default=50)
    parser.add_argument("--global-seed", type=int, default=42)
    parser.add_argument("--min-cluster-purity", type=float, default=0.9)
    parser.add_argument("--min-rf-accuracy", type=float, default=0.9)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    np.random.seed(args.global_seed)

    MIN_SAMPLES = args.min_samples
    MIN_CLUSTER_PURITY = args.min_cluster_purity
    MIN_RF_ACCURACY = args.min_rf_accuracy

    OUTPUT_DIR = Path(args.out_dir)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"[gtex] Output directory: {OUTPUT_DIR}", flush=True)

    # Load data
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

    # transpose to get samples as rows, LVs as columns
    lv_matrix = lv_data.T

    # filter metadata to only include samples in LV matrix
    meta_filtered = gtex_meta[gtex_meta["SAMPID"].isin(lv_matrix.index)].copy()
    meta_filtered = meta_filtered.set_index("SAMPID")
    meta_filtered = meta_filtered.loc[lv_matrix.index]

    # filter to tissues with >= MIN_SAMPLES samples
    tissue_counts = meta_filtered["SMTSD"].value_counts()
    valid_tissues = tissue_counts[tissue_counts >= MIN_SAMPLES].index
    mask = meta_filtered["SMTSD"].isin(valid_tissues)
    meta_filtered = meta_filtered[mask]
    lv_matrix = lv_matrix.loc[meta_filtered.index]

    tissue_labels = meta_filtered["SMTSD"]

    print(f"[gtex] Number of samples: {len(lv_matrix)}", flush=True)
    print(f"[gtex] Number of LVs: {lv_matrix.shape[1]}", flush=True)
    print(f"[gtex] Number of unique tissues (>= {MIN_SAMPLES} samples): {tissue_labels.nunique()}", flush=True)

    # Load K-means model
    with open(args.kmeans_model, "rb") as f:
        model_bundle = pickle.load(f)

    km_model = model_bundle["model"]
    km_scaler = model_bundle["scaler"]
    best_k = model_bundle["best_k"]
    best_approach = model_bundle["best_approach"]
    saved_ari = model_bundle["ari"]

    print(f"[gtex] K-means model: k={best_k}, approach={best_approach}, ARI={saved_ari:.4f}", flush=True)

    # apply the same preprocessing and get cluster assignments
    if best_approach == "scaled":
        X_km = km_scaler.transform(lv_matrix.astype(np.float32))
    else:
        X_km = lv_matrix.values.astype(np.float32)

    cluster_labels = km_model.predict(X_km)

    # verify ARI matches
    reproduced_ari = adjusted_rand_score(tissue_labels.values, cluster_labels)
    print(f"[gtex] Reproduced ARI: {reproduced_ari:.4f} (saved: {saved_ari:.4f})", flush=True)

    # map clusters to tissues (majority vote)
    cluster_tissue_map = {}
    cluster_tissue_purity = {}

    for c in range(best_k):
        cmask = cluster_labels == c
        cluster_tissues = tissue_labels[cmask]
        dominant_tissue = cluster_tissues.value_counts().index[0]
        purity = cluster_tissues.value_counts().iloc[0] / len(cluster_tissues)
        cluster_tissue_map[c] = dominant_tissue
        cluster_tissue_purity[c] = purity

    # group clusters by dominant tissue
    tissue_to_clusters = defaultdict(list)
    for c, tissue in cluster_tissue_map.items():
        tissue_to_clusters[tissue].append(c)

    cluster_summary = pd.DataFrame({
        "Cluster": range(best_k),
        "Dominant_Tissue": [cluster_tissue_map[c] for c in range(best_k)],
        "Purity": [cluster_tissue_purity[c] for c in range(best_k)],
        "Size": [np.sum(cluster_labels == c) for c in range(best_k)],
    }).sort_values("Dominant_Tissue")
    print(f"[gtex] Cluster summary ({best_k} clusters):", flush=True)
    print(cluster_summary, flush=True)

    # binary tissue labels derived from k-means (positive = sample in a cluster dominated by that tissue)
    km_tissue_labels = pd.Series(
        [cluster_tissue_map[c] for c in cluster_labels],
        index=lv_matrix.index,
    )
    print(f"[gtex] K-means-derived tissue label distribution:", flush=True)
    print(km_tissue_labels.value_counts(), flush=True)

    # Feature importance analysis (K-means binary + SHAP)
    tissues = sorted(tissue_to_clusters.keys())
    tissues = [
        t for t in tissues
        if np.mean([cluster_tissue_purity[c] for c in tissue_to_clusters[t]]) >= MIN_CLUSTER_PURITY
    ]
    print(f"[gtex] Tissues with purity >= {MIN_CLUSTER_PURITY}: {tissues}", flush=True)

    def extra_accuracy_fields(tissue: str) -> dict:
        clusters = tissue_to_clusters[tissue]
        return {
            "N_Clusters": len(clusters),
            "Mean_Purity": np.mean([cluster_tissue_purity[c] for c in clusters]),
        }

    run_tissue_rf_shap(
        lv_matrix,
        tissues,
        label_fn=lambda tissue: (km_tissue_labels == tissue),
        output_dir=OUTPUT_DIR,
        global_seed=args.global_seed,
        min_rf_accuracy=MIN_RF_ACCURACY,
        param_distributions=PARAM_DISTRIBUTIONS,
        extra_accuracy_fields_fn=extra_accuracy_fields,
    )


if __name__ == "__main__":
    main()
