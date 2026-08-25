#!/usr/bin/env python3
# GTEx model building with PCA, NMF and ICA.

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA, NMF, FastICA


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--df-gtex-fbm-filt", required=True)
    parser.add_argument("--k", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument(
        "--method",
        choices=("all", "PCA", "NMF", "ICA"),
        default="all",
        help="Fit one method or all three (default: all, the production behavior). "
        "The timing benchmark fits one per process so that ICA and NMF are not "
        "charged the preceding PCA fit and CSV parse.",
    )
    parser.add_argument(
        "--flat-output",
        action="store_true",
        help="Write a single selected method straight into --out-dir instead of a "
        "method subdirectory.",
    )
    args = parser.parse_args()

    gtex_data = pd.read_csv(args.df_gtex_fbm_filt, index_col=0).astype(np.float32)

    K = pd.read_csv(args.k)
    n_components = int(K["CLAMP_K_gtex"].iloc[0])
    print(f"[gtex] Number of components: {n_components}", flush=True)

    root = Path(args.out_dir)

    def out_for(name: str) -> Path:
        directory = root if args.flat_output else root / name
        directory.mkdir(parents=True, exist_ok=True)
        return directory

    if args.method in ("all", "PCA"):
        pca_dir = out_for("PCA")
        X = gtex_data.T

        pca = PCA(n_components=n_components, svd_solver="auto", random_state=args.seed)
        W = pca.fit_transform(X)
        H = pca.components_

        pc_names = [f"PC{i + 1}" for i in range(W.shape[1])]

        gtex_pca_scores = pd.DataFrame(W, index=gtex_data.columns, columns=pc_names)
        gtex_pca_loadings = pd.DataFrame(H, index=pc_names, columns=gtex_data.index)
        gtex_pca_B = gtex_pca_scores.T
        gtex_pca_B.index.name = "PC"

        gtex_pca_B.to_pickle(pca_dir / "gtex_pca_B.pkl")
        gtex_pca_scores.to_pickle(pca_dir / "gtex_pca_scores.pkl")
        gtex_pca_loadings.to_pickle(pca_dir / "gtex_pca_loadings.pkl")
        print(f"[gtex] PCA saved -> {pca_dir}", flush=True)

    if args.method in ("all", "ICA"):
        ica_dir = out_for("ICA")
        X = gtex_data.T

        ica = FastICA(n_components=n_components, random_state=args.seed, max_iter=2000)
        W = ica.fit_transform(X)
        H = ica.mixing_.T

        ic_names = [f"IC{i + 1}" for i in range(W.shape[1])]

        gtex_ica_scores = pd.DataFrame(W, index=gtex_data.columns, columns=ic_names)
        gtex_ica_loadings = pd.DataFrame(H, index=ic_names, columns=gtex_data.index)
        gtex_ica_B = gtex_ica_scores.T
        gtex_ica_B.index.name = "IC"

        gtex_ica_B.to_pickle(ica_dir / "gtex_ica_B.pkl")
        gtex_ica_scores.to_pickle(ica_dir / "gtex_ica_scores.pkl")
        gtex_ica_loadings.to_pickle(ica_dir / "gtex_ica_loadings.pkl")
        print(f"[gtex] ICA saved -> {ica_dir}", flush=True)

    if args.method in ("all", "NMF"):
        nmf_dir = out_for("NMF")
        gene_min = gtex_data.min(axis=1)
        gtex_data_nmf = gtex_data.sub(gene_min, axis=0)

        X = gtex_data_nmf.T

        nmf = NMF(n_components=n_components, init="nndsvd", random_state=args.seed, max_iter=1000)
        W = nmf.fit_transform(X)
        H = nmf.components_

        lv_names = [f"LV{i + 1}" for i in range(W.shape[1])]

        gtex_nmf_scores = pd.DataFrame(W, index=gtex_data.columns, columns=lv_names)
        gtex_nmf_loadings = pd.DataFrame(H, index=lv_names, columns=gtex_data.index)
        gtex_nmf_B = gtex_nmf_scores.T
        gtex_nmf_B.index.name = "LV"

        gtex_nmf_B.to_pickle(nmf_dir / "gtex_nmf_B.pkl")
        gtex_nmf_scores.to_pickle(nmf_dir / "gtex_nmf_scores.pkl")
        gtex_nmf_loadings.to_pickle(nmf_dir / "gtex_nmf_loadings.pkl")
        print(f"[gtex] NMF saved -> {nmf_dir}", flush=True)


if __name__ == "__main__":
    main()
