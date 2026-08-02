#!/usr/bin/env python3
"""Derive the cohort sample manifest for a raw-built pseudobulk dataset."""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from common import obs_values, read_yaml, repo_path, transform_ids


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    config = read_yaml(args.config)
    cfg = config["datasets"][args.dataset]
    if cfg["kind"] != "h5ad":
        raise ValueError(f"cohort manifest generation only supports h5ad datasets, got {cfg['kind']!r}")

    import h5py

    # Map each cell to its sample, then blank out cells excluded by config
    with h5py.File(repo_path(cfg["raw"]), "r") as h5:
        eligible_sample = transform_ids(obs_values(h5, cfg["sample_col"]), cfg.get("sample_transform"))
        for column, excluded_values in cfg.get("exclude", {}).items():
            excluded = np.isin(obs_values(h5, column), np.asarray(excluded_values).astype(str))
            eligible_sample[excluded] = ""

    # Count cells per sample and drop samples below the min_cells threshold
    counts = pd.Series(eligible_sample).value_counts()
    counts = counts[counts.index != ""]
    min_cells = cfg.get("min_cells")
    if min_cells is not None:
        counts = counts[counts >= min_cells]
    manifest = counts.rename_axis("sample").reset_index(name="nCells").sort_values("sample")

    repo_path(args.output).parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(repo_path(args.output), index=False)
    print(f"[{args.dataset}] cohort manifest: {len(manifest)} samples")


if __name__ == "__main__":
    main()
