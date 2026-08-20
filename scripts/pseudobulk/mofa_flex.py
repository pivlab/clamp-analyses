#!/usr/bin/env python3
"""MOFA-FLEX GO:BP pseudobulk models"""

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import anndata as ad
import mofaflex as mfl
from mofaflex._core.feature_sets import FeatureSets
import pandas as pd
import torch


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--norm", required=True)
    parser.add_argument("--k", required=True)
    parser.add_argument("--gmt", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--max-epochs", type=int, default=200)
    parser.add_argument("--seed", type=int, default=123)
    parser.add_argument(
        "--threads",
        type=int,
        default=None,
        help="Optionally cap PyTorch intra-op and inter-op thread pools.",
    )
    args = parser.parse_args()
    if args.threads is not None:
        torch.set_num_threads(args.threads)
        torch.set_num_interop_threads(args.threads)
    norm = pd.read_csv(args.norm, index_col=0).astype("float32")
    k = int(pd.read_csv(args.k)["k"].iloc[0])
    genes, samples = norm.index.astype(str).tolist(), norm.columns.astype(str).tolist()
    collection = FeatureSets.from_gmt(args.gmt, name="GO_Biological_Process_pinned")
    collection = collection.filter(genes, min_fraction=0.4, min_count=40, max_count=200)
    collection = collection.merge_similar(metric="jaccard", similarity_threshold=0.8, iteratively=True)
    data = ad.AnnData(X=norm.T.values, obs=pd.DataFrame(index=samples), var=pd.DataFrame(index=genes))
    data.varm["annotations"] = collection.to_mask(genes).T

    model = mfl.MOFAFLEX(
        {"group_1": {"view_1": data}},
        mfl.DataOptions(scale_per_group=False, plot_data_overview=False, annotations_varm_key="annotations"),
        mfl.ModelOptions(n_factors=k, weight_prior="Horseshoe", likelihoods="Normal"),
        mfl.TrainingOptions(seed=args.seed, max_epochs=args.max_epochs, save_path=False, device="cpu"),
    )
    factors = model.get_factors()["group_1"]
    weights = model.get_weights()["view_1"]

    factors = factors.T
    factors.columns = samples
    factors.index = [f"LV{i + 1}" for i in range(len(factors))]
    weights = weights.T
    weights.index = genes
    weights.columns = factors.index
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    factors.to_csv(out / "B.csv")
    weights.to_csv(out / "Z.csv")
    with open(out / "model.pkl", "wb") as handle:
        pickle.dump(model, handle)


if __name__ == "__main__":
    main()
