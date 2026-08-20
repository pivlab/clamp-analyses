#!/usr/bin/env python3
"""Fit PCA, NMF, and ICA to one preprocessed pseudobulk matrix."""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import FastICA, NMF, PCA


def save_model(root: Path, method: str, scores: np.ndarray, loadings: np.ndarray,
               sample_names: list[str], gene_names: list[str], model: object,
               prefix: str, filename: str, flat_output: bool = False) -> None:
    out = root if flat_output else root / method
    out.mkdir(parents=True, exist_ok=True)
    names = [f"{prefix}{i + 1}" for i in range(scores.shape[1])]
    # LVs x samples / genes x LVs, matching CLAMP's B/Z convention
    pd.DataFrame(scores.T, index=names, columns=sample_names).to_csv(out / "B.csv")
    pd.DataFrame(loadings.T, index=gene_names, columns=names).to_csv(out / "Z.csv")
    with open(out / filename, "wb") as handle:
        pickle.dump(model, handle)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--norm", required=True)
    parser.add_argument("--k", required=True)
    parser.add_argument("--models-dir", required=True)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument(
        "--method",
        choices=("all", "PCA", "NMF", "ICA"),
        default="all",
        help="Fit one method or all three (default: all, preserving production behavior).",
    )
    parser.add_argument(
        "--flat-output",
        action="store_true",
        help="Write a selected method directly into --models-dir instead of a method subdirectory.",
    )
    args = parser.parse_args()
    if args.flat_output and args.method == "all":
        parser.error("--flat-output requires one explicit --method (PCA, NMF, or ICA)")

    norm = pd.read_csv(args.norm, index_col=0).astype(np.float32)
    k = int(pd.read_csv(args.k)["k"].iloc[0])
    samples, genes = norm.columns.astype(str).tolist(), norm.index.astype(str).tolist()
    x = norm.T.to_numpy()
    root = Path(args.models_dir)

    if args.method in ("all", "PCA"):
        pca = PCA(n_components=k, random_state=args.seed)
        pca_scores = pca.fit_transform(x)
        save_model(root, "PCA", pca_scores, pca.components_, samples, genes, pca,
                   "PC", "pca_model.pkl", args.flat_output)

    if args.method in ("all", "ICA"):
        ica = FastICA(n_components=k, random_state=args.seed, max_iter=2000, tol=0.01)
        ica_scores = ica.fit_transform(x)
        save_model(root, "ICA", ica_scores, ica.mixing_.T, samples, genes, ica,
                   "IC", "ica_model.pkl", args.flat_output)

    if args.method in ("all", "NMF"):
        # NMF needs non-negative input, so shift each gene's z-scores by its row min
        nonnegative = norm.sub(norm.min(axis=1), axis=0).T.to_numpy()
        nmf = NMF(n_components=k, init="nndsvd", random_state=args.seed, max_iter=1000)
        nmf_scores = nmf.fit_transform(nonnegative)
        save_model(root, "NMF", nmf_scores, nmf.components_, samples, genes, nmf,
                   "LV", "nmf_model.pkl", args.flat_output)


if __name__ == "__main__":
    main()
