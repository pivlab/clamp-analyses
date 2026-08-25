#!/usr/bin/env python3
# MOFA-FLEX model on GTEx data with Biological Process (BP) pathway priors.

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import anndata as ad
import mofaflex as mfl
import numpy as np
import pandas as pd
import torch
from mofaflex._core.feature_sets import FeatureSets


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--df-gtex-fbm-filt", required=True)
    parser.add_argument("--k", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--gmt", required=True)
    parser.add_argument("--min-fraction", type=float, default=0.4)
    parser.add_argument("--min-count", type=int, default=40)
    parser.add_argument("--max-count", type=int, default=200)
    parser.add_argument("--similarity-threshold", type=float, default=0.8)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument("--max-epochs", type=int, default=1000)
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Cap the PyTorch intra-op and inter-op thread pools.  torch otherwise "
        "defaults to one thread per core, which would give MOFA-FLEX several times "
        "the parallelism of every other method in the timing benchmark.",
    )
    args = parser.parse_args()

    if args.threads is not None:
        torch.set_num_threads(args.threads)
        torch.set_num_interop_threads(args.threads)

    gtex_data = pd.read_csv(args.df_gtex_fbm_filt, index_col=0).astype(np.float32)
    print(f"[gtex] GTEx data shape: {gtex_data.shape}", flush=True)
    print(f"[gtex]   Samples: {gtex_data.shape[0]}", flush=True)
    print(f"[gtex]   Genes: {gtex_data.shape[1]}", flush=True)

    K = pd.read_csv(args.k)
    n_components = int(K["CLAMP_K_gtex"].iloc[0])
    print(f"[gtex] Number of components: {n_components}", flush=True)

    gene_list = gtex_data.index.tolist()

    bp_collection = FeatureSets.from_gmt(args.gmt, name="GO_Biological_Process_pinned")

    bp_collection = bp_collection.filter(
        gene_list,
        min_fraction=args.min_fraction,
        min_count=args.min_count,
        max_count=args.max_count,
    )

    bp_collection = bp_collection.merge_similar(
        metric="jaccard",
        similarity_threshold=args.similarity_threshold,
        iteratively=True,
    )
    print(f"[gtex] Filtered pathways: {len(bp_collection)}", flush=True)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    adata = ad.AnnData(
        X=gtex_data.T.values,
        obs=pd.DataFrame(index=gtex_data.columns),
        var=pd.DataFrame(index=gtex_data.index),
    )
    print(f"[gtex] AnnData shape: {adata.shape}", flush=True)
    print(f"[gtex]   Samples (obs): {adata.n_obs}", flush=True)
    print(f"[gtex]   Genes (var): {adata.n_vars}", flush=True)

    adata.varm["annotations"] = bp_collection.to_mask(gene_list).T
    print(f"[gtex] Annotations shape: {adata.varm['annotations'].shape}", flush=True)
    print(f"[gtex] Number of BP pathways: {adata.varm['annotations'].shape[1]}", flush=True)

    data_options = mfl.DataOptions(
        scale_per_group=False,
        plot_data_overview=False,
        annotations_varm_key="annotations",
    )
    model_options = mfl.ModelOptions(
        n_factors=n_components,
        weight_prior="Horseshoe",
        likelihoods="Normal",
    )
    training_options = mfl.TrainingOptions(
        seed=args.seed,
        max_epochs=args.max_epochs,
        save_path=False,
        device="cpu",
    )

    model = mfl.MOFAFLEX(
        {"group_1": {"view_1": adata}},
        data_options,
        model_options,
        training_options,
    )

    factors = model.get_factors()["group_1"]
    B_matrix = factors.T

    weights = model.get_weights()["view_1"]

    B_matrix.to_csv(out_dir / "B_matrix.csv")

    weights.to_csv(out_dir / "Z_matrix.csv")

    with open(out_dir / "model.pkl", "wb") as f:
        pickle.dump(model, f)
    print(f"[gtex] Saved model to {out_dir / 'model.pkl'}", flush=True)


if __name__ == "__main__":
    main()
