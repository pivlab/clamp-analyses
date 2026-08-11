#!/usr/bin/env python3
"""Build an independent expression-derived UMAP for sampled retained cells."""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


REPO_ROOT = Path(__file__).resolve().parents[2]


def path(value: str) -> Path:
    candidate = Path(value)
    return candidate if candidate.is_absolute() else REPO_ROOT / candidate


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--cells", required=True)
    parser.add_argument("--genes", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--target-sum", type=float, default=10_000)
    parser.add_argument("--n-hvg", type=int, default=2_000)
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--min-dist", type=float, default=0.30)
    parser.add_argument("--seed", type=int, default=123)
    args = parser.parse_args()

    counts = sp.load_npz(path(args.counts)).tocsr().astype(np.float32)
    cells = pd.read_csv(path(args.cells), dtype={"cell_id": str, "sample": str})
    genes = pd.read_csv(path(args.genes), dtype=str)["gene"].astype(str)
    if counts.shape != (len(cells), len(genes)):
        raise ValueError(f"{args.dataset}: sampled matrix and metadata dimensions differ")
    if (
        cells["cell_index"].duplicated().any()
        or cells["cell_id"].duplicated().any()
        or genes.duplicated().any()
    ):
        raise ValueError(f"{args.dataset}: sampled cells or genes are not unique")
    if (
        not np.isfinite(counts.data).all()
        or np.any(counts.data < 0)
        or not np.allclose(counts.data, np.round(counts.data))
    ):
        raise ValueError(f"{args.dataset}: sampled matrix contains invalid values")

    adata = ad.AnnData(X=counts, obs=cells.set_index("cell_id", drop=False))
    adata.var_names = genes.to_numpy()
    sc.pp.normalize_total(adata, target_sum=args.target_sum)
    sc.pp.log1p(adata)
    n_hvg = min(args.n_hvg, adata.n_vars)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, flavor="seurat")
    if int(adata.var["highly_variable"].sum()) < 2:
        raise ValueError(f"{args.dataset}: fewer than two highly variable genes")
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    n_pcs = min(args.n_pcs, adata.n_obs - 1, adata.n_vars - 1)
    if n_pcs < 2:
        raise ValueError(f"{args.dataset}: sampled expression matrix cannot support PCA")
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver="arpack", random_state=args.seed)
    n_neighbors = min(args.n_neighbors, adata.n_obs - 1)
    sc.pp.neighbors(
        adata, n_neighbors=n_neighbors, n_pcs=n_pcs,
        method="umap", random_state=args.seed,
    )
    sc.tl.umap(adata, min_dist=args.min_dist, random_state=args.seed)
    coords = np.asarray(adata.obsm["X_umap"], dtype=np.float64)
    if coords.shape != (len(cells), 2) or not np.isfinite(coords).all():
        raise ValueError(f"{args.dataset}: UMAP output is incomplete or non-finite")

    output = cells.copy()
    output["dataset"] = args.dataset
    output["umap1"] = coords[:, 0]
    output["umap2"] = coords[:, 1]
    output_path = path(args.output)
    summary_path = path(args.summary)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, index=False)
    pd.DataFrame(
        [
            {
                "dataset": args.dataset,
                "n_cells": adata.n_obs,
                "n_genes_input": len(genes),
                "n_hvg": adata.n_vars,
                "n_pcs": n_pcs,
                "n_neighbors": n_neighbors,
                "target_sum": args.target_sum,
                "min_dist": args.min_dist,
                "seed": args.seed,
                "embedding_source": "sampled_raw_expression_normalize_hvg_pca_neighbors_umap",
                "all_coordinates_finite": True,
            }
        ]
    ).to_csv(summary_path, index=False)
    print(f"[{args.dataset}] wrote expression UMAP for {adata.n_obs:,} cells")


if __name__ == "__main__":
    main()
