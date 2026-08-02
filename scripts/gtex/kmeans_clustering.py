#!/usr/bin/env python3
"""GTEx k-means clustering of gene expression data using different methods.

GPU-accelerated (cuML) k-means ensemble over sample embeddings from multiple
model-building methods plus raw RNA-Seq, scored against true tissue labels
(ARI). Fitting half of nbs/03_model_biology/01_gtex/00_kmeans_clustering.ipynb;
the notebook itself only reads this script's cached results and plots them.
"""

from __future__ import annotations

import argparse
import csv
import pickle
import sys
from pathlib import Path

import cupy as cp
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from cuml.cluster import KMeans as cuKMeans
from sklearn.cluster import KMeans as skKMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import extract_B_matrix  # noqa: E402

readRDS = ro.r["readRDS"]


def load_rds_B_matrix(path: str) -> pd.DataFrame:
    return extract_B_matrix(readRDS(str(path)))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtex-meta", required=True)
    parser.add_argument("--df-gtex-fbm-filt", required=True)
    parser.add_argument("--clamp-base-rds", required=True)
    parser.add_argument("--clamp-full-rds", required=True)
    parser.add_argument("--plier-rds", required=True)
    parser.add_argument("--gss-b", required=True)
    parser.add_argument("--flashier-b", required=True)
    parser.add_argument("--mofa-flex-prior-b", required=True)
    parser.add_argument("--pca-b", required=True)
    parser.add_argument("--ica-b", required=True)
    parser.add_argument("--nmf-b", required=True)
    parser.add_argument("--min-samples", type=int, default=50)
    parser.add_argument("--base-seed", type=int, default=10000)
    parser.add_argument("--n-reps-per-k", type=int, default=50)
    parser.add_argument(
        "--gene-fractions", type=float, nargs="+",
        default=[1.0, 0.75, 0.50, 0.25, 0.10, 0.05, 0.01],
    )
    parser.add_argument("--cache-dir", required=True)
    parser.add_argument("--model-out-dir", required=True)
    args = parser.parse_args()

    MIN_SAMPLES = args.min_samples
    BASE_SEED = args.base_seed
    N_REPS_PER_K = args.n_reps_per_k
    GENE_FRACTIONS = args.gene_fractions

    KMEANS_CACHE_DIR = Path(args.cache_dir)
    KMEANS_CACHE_DIR.mkdir(parents=True, exist_ok=True)
    model_dir = Path(args.model_out_dir)
    model_dir.mkdir(parents=True, exist_ok=True)

    # Input
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

    total_tissues = len(gtex_meta["SMTSD"].value_counts())

    gtex_data = pd.read_csv(args.df_gtex_fbm_filt, index_col=0)

    meta = gtex_meta.copy()
    if meta.index.name != "SAMPID":
        meta = meta.set_index("SAMPID")

    sample_ids_expr = pd.Index(gtex_data.columns.astype(str).str.strip(), name="SAMPID")
    sample_ids_meta = meta.index.astype(str).str.strip()
    common_ids = sample_ids_expr.intersection(sample_ids_meta)
    if len(common_ids) == 0:
        raise ValueError("No overlapping sample IDs between expression and metadata.")

    gtex_data_common = gtex_data.loc[:, common_ids]
    meta_common = meta.loc[common_ids]

    tissue_counts = meta_common["SMTSD"].value_counts()
    valid_tissues = tissue_counts[tissue_counts >= MIN_SAMPLES].index
    mask = meta_common["SMTSD"].isin(valid_tissues)
    meta_common = meta_common[mask]
    common_ids = meta_common.index
    gtex_data_common = gtex_data_common.loc[:, common_ids]

    total_tissues = len(valid_tissues)
    print(f"[gtex] Tissues with >= {MIN_SAMPLES} samples: {total_tissues}", flush=True)

    # GPU kmeans ensemble
    k_values = [total_tissues]

    def load_or_run_cached(cache_key, compute_fn):
        cache_path = KMEANS_CACHE_DIR / f"{cache_key}.pkl"
        if cache_path.exists():
            print(f"[gtex] Loading cached results from {cache_path}", flush=True)
            with open(cache_path, "rb") as f:
                return pickle.load(f)

        results = compute_fn()

        with open(cache_path, "wb") as f:
            pickle.dump(results, f)
        print(f"[gtex] Saved results to {cache_path}", flush=True)
        return results

    def run_kmeans_compare_scaling(B_t, y_true, method_name, emb_k_values=None):
        if emb_k_values is None:
            emb_k_values = k_values

        X = StandardScaler().fit_transform(B_t)
        X_gpu = cp.asarray(X.astype("float32"))

        records = []
        for k in emb_k_values:
            for rep in range(N_REPS_PER_K):
                km = cuKMeans(
                    n_clusters=k,
                    init="k-means++",
                    random_state=BASE_SEED + k * 1000 + rep,
                    max_iter=300,
                    tol=1e-4,
                    verbose=0,
                )
                labels_gpu = km.fit_predict(X_gpu)
                labels = cp.asnumpy(labels_gpu)
                ari = adjusted_rand_score(y_true, labels)

                records.append({
                    "k": k,
                    "rep": rep,
                    "seed": BASE_SEED + k * 1000 + rep,
                    "ari": ari,
                })

        del X_gpu
        cp.get_default_memory_pool().free_all_blocks()

        df = pd.DataFrame.from_records(records)
        mean_ari = df.groupby("k")["ari"].mean()
        best_k = mean_ari.idxmax()
        best_k_df = df[df["k"] == best_k]

        print(f"[gtex] {method_name}: Best k: {best_k}, Mean ARI: {mean_ari[best_k]:.4f}", flush=True)

        return {
            "best_results_df": df,
            "best_k_results_df": best_k_df,
            "best_k": best_k,
            "scaled_ari": mean_ari[best_k],
        }

    def prepare_embedding_for_kmeans(B, meta):
        B = B.copy()
        B.columns = B.columns.astype(str).str.strip()

        common_ids = B.columns.intersection(meta.index)
        if len(common_ids) == 0:
            raise ValueError("No overlapping sample IDs between B and metadata.")

        meta_f = meta.loc[common_ids]
        tissue_counts = meta_f["SMTSD"].value_counts()
        valid_tissues = tissue_counts[tissue_counts >= MIN_SAMPLES].index
        mask = meta_f["SMTSD"].isin(valid_tissues)
        meta_f = meta_f[mask]
        common_ids = meta_f.index

        B_t = B.loc[:, common_ids].T.astype(np.float32)
        y_true = meta_f["SMTSD"].astype(str).to_numpy()
        n_tissues = len(valid_tissues)
        return B_t, y_true, n_tissues

    def run_embedding_kmeans(B, meta, method_name, cache_key):
        B_t, y_true, n_tissues = prepare_embedding_for_kmeans(B, meta)
        emb_k_values = [n_tissues]
        results = load_or_run_cached(
            cache_key,
            lambda: run_kmeans_compare_scaling(B_t, y_true, method_name, emb_k_values),
        )
        print(f"[gtex] {method_name}: Best k: {results['best_k']}, Mean ARI: {results['scaled_ari']:.4f}", flush=True)
        return results

    # RNA-Seq: gene-selection-by-variance sweep
    def compute_rnaseq_topvar_results():
        B = gtex_data_common.copy()
        B.columns = B.columns.astype(str).str.strip()

        common_ids = B.columns.intersection(meta_common.index)
        B_common = B.loc[:, common_ids]

        n_genes = len(B_common.index)
        y_true = meta_common.loc[common_ids, "SMTSD"].astype(str).to_numpy()
        rnaseq_k_values = [total_tissues]

        gene_variances = B_common.var(axis=1)
        genes_by_variance = gene_variances.sort_values(ascending=False).index.tolist()

        print(f"[gtex] Total genes: {n_genes}", flush=True)
        print(f"[gtex] Total samples: {len(common_ids)}", flush=True)
        print(f"[gtex] k = {total_tissues} tissues", flush=True)

        results_by_fraction = {}

        for frac in GENE_FRACTIONS:
            pct = int(frac * 100)
            run_label = f"rnaseq_{pct}pct_topvar"
            print(f"[gtex] Running: {run_label}", flush=True)

            n_select = int(n_genes * frac)
            selected_genes = genes_by_variance[:n_select]

            B_subset = B_common.loc[selected_genes, :]
            B_t = B_subset.T.astype(np.float32)

            results_by_fraction[run_label] = run_kmeans_compare_scaling(B_t, y_true, run_label, rnaseq_k_values)

        for label, results in results_by_fraction.items():
            print(f"[gtex]   {label}: Best k = {results['best_k']}, ARI = {results['scaled_ari']:.4f}", flush=True)

        return results_by_fraction

    rnaseq_topvar_results = load_or_run_cached(
        "gtex_rnaseq_topvar_results",
        compute_rnaseq_topvar_results,
    )

    # CLAMP models
    CLAMP_RDS_FILES = {
        "CLAMPbase": {"model_path": args.clamp_base_rds},
        "CLAMPfull": {"model_path": args.clamp_full_rds},
    }

    clamp_embeddings = {}
    for model_name, spec in CLAMP_RDS_FILES.items():
        clamp_embeddings[model_name] = load_rds_B_matrix(spec["model_path"])
        print(f"[gtex] {model_name}: {clamp_embeddings[model_name].shape}", flush=True)

    gtex_CLAMP_base_results = run_embedding_kmeans(
        clamp_embeddings["CLAMPbase"], meta, "CLAMP_base", "gtex_CLAMP_base_kmeans",
    )
    gtex_CLAMPfull_results = run_embedding_kmeans(
        clamp_embeddings["CLAMPfull"], meta, "CLAMPfull", "gtex_CLAMPfull_kmeans",
    )

    # Final CLAMPfull model: refit + save a portable sklearn-compatible bundle
    B_clamp_full = clamp_embeddings["CLAMPfull"].copy()
    B_clamp_full.columns = B_clamp_full.columns.astype(str).str.strip()
    common_ids_clamp = B_clamp_full.columns.intersection(meta_common.index)
    B_t_clamp = B_clamp_full.loc[:, common_ids_clamp].T.astype(np.float32)

    best_k = gtex_CLAMPfull_results["best_k"]

    scaler = StandardScaler()
    X_clamp = scaler.fit_transform(B_t_clamp)

    best_row = gtex_CLAMPfull_results["best_k_results_df"].loc[
        gtex_CLAMPfull_results["best_k_results_df"]["ari"].idxmax()
    ]
    best_seed = int(best_row["seed"])
    final_ari = float(best_row["ari"])

    X_gpu = cp.asarray(X_clamp.astype("float32"))
    final_km = cuKMeans(
        n_clusters=best_k,
        init="k-means++",
        random_state=best_seed,
        max_iter=300,
        tol=1e-4,
        verbose=0,
    )
    final_km.fit(X_gpu)

    # Convert to sklearn so downstream notebooks don't need cuml
    centers_np = cp.asnumpy(final_km.cluster_centers_)
    sklearn_km = skKMeans(n_clusters=best_k, init=centers_np, n_init=1, max_iter=1)
    sklearn_km.fit(centers_np)
    sklearn_km.cluster_centers_ = centers_np

    del X_gpu
    cp.get_default_memory_pool().free_all_blocks()

    model_bundle = {
        "model": sklearn_km,
        "scaler": scaler,
        "best_approach": "scaled",
        "best_k": best_k,
        "best_seed": best_seed,
        "ari": final_ari,
        "feature_names": list(B_t_clamp.columns),
        "sample_ids": list(B_t_clamp.index),
    }

    save_path = model_dir / "gtex_CLAMPfull_kmeans_model.pkl"
    with open(save_path, "wb") as f:
        pickle.dump(model_bundle, f)

    print(f"[gtex] Model saved to: {save_path}", flush=True)
    print(f"[gtex]   best_k={best_k}, seed={best_seed}, ARI={final_ari:.4f}", flush=True)

    # PLIER
    gtex_PLIER_B = load_rds_B_matrix(args.plier_rds)
    run_embedding_kmeans(gtex_PLIER_B, meta, "PLIER", "gtex_PLIER_kmeans")

    # GenomicSuperSignature
    gtex_GenomicSuperSignature_B = pd.read_csv(args.gss_b, index_col=0)
    run_embedding_kmeans(
        gtex_GenomicSuperSignature_B, meta, "GenomicSuperSignature", "gtex_GenomicSuperSignature_kmeans",
    )

    # Flashier
    gtex_flashier_B = pd.read_csv(args.flashier_b, index_col=0)
    run_embedding_kmeans(gtex_flashier_B, meta, "Flashier", "gtex_flashier_kmeans")

    # MOFA FLEX BP
    gtex_mofabc_B = pd.read_csv(args.mofa_flex_prior_b, index_col=0)
    run_embedding_kmeans(gtex_mofabc_B, meta, "MOFA-FLEX_BP", "gtex_mofabc_kmeans")

    # PCA
    gtex_pca_B = pd.read_pickle(args.pca_b)
    run_embedding_kmeans(gtex_pca_B, meta, "PCA", "gtex_pca_kmeans")

    # ICA
    gtex_ica_B = pd.read_pickle(args.ica_b)
    run_embedding_kmeans(gtex_ica_B, meta, "ICA", "gtex_ica_kmeans")

    # NMF
    gtex_nmf_B = pd.read_pickle(args.nmf_b)
    run_embedding_kmeans(gtex_nmf_B, meta, "NMF", "gtex_nmf_kmeans")

    print("[gtex] kmeans_clustering done", flush=True)


if __name__ == "__main__":
    main()
