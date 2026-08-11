#!/usr/bin/env python3
"""End-to-end integrity checks for the donor-bulk benchmark."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]


def resolve(value: str | Path) -> Path:
    value = Path(value)
    return value if value.is_absolute() else REPO_ROOT / value


def matrix(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    return frame.astype(float)


def finite_projection(path: Path, chunk: int = 50_000) -> tuple[int, int]:
    with h5py.File(path, "r") as h5:
        scores = h5["scores"]
        for start in range(0, scores.shape[0], chunk):
            if not np.isfinite(scores[start : start + chunk]).all():
                raise ValueError(f"non-finite projected scores: {path}")
        if h5.attrs.get("duplicate_gene_policy", "") != "sum":
            raise ValueError(f"projection does not record donor-sum duplicate semantics: {path}")
        return scores.shape


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--production-root", required=True)
    parser.add_argument("--analysis-dir", required=True)
    parser.add_argument("--output-summary", required=True)
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--datasets", nargs="+", required=True)
    parser.add_argument("--sample-target", type=int, default=60_000)
    args = parser.parse_args()
    production = resolve(args.production_root)
    analysis = resolve(args.analysis_dir)
    rows = []
    for dataset in args.datasets:
        root = production / dataset
        aggregation = pd.read_csv(root / "bulk/aggregation_summary.csv").iloc[0]
        if aggregation.raw_validation_scope != "exhaustive_all_nonzeros":
            raise ValueError(f"{dataset}: raw input was not exhaustively validated")
        donor = matrix(root / "bulk/bulk_counts.csv")
        control = matrix(root / "bulk/bulk_mean_cell_cpm.csv")
        truth = matrix(root / "bulk/truthFrac_v0.csv")
        info = pd.read_csv(root / "bulk/patient_info.csv", dtype={"sample": str})
        donors = info["sample"].astype(str).tolist()
        if donor.index.duplicated().any() or donor.columns.duplicated().any():
            raise ValueError(f"{dataset}: donor matrix has duplicate identifiers")
        if not donor.index.equals(control.index):
            raise ValueError(f"{dataset}: donor/control genes differ")
        if list(donor.columns) != donors or list(control.columns) != donors or list(truth.index) != donors:
            raise ValueError(f"{dataset}: donor alignment failure")
        values = donor.to_numpy()
        if not np.isfinite(values).all() or np.any(values < 0) or not np.allclose(values, np.round(values)):
            raise ValueError(f"{dataset}: invalid donor raw sums")
        if not np.allclose(control.sum(axis=0), 1e6, atol=1e-3):
            raise ValueError(f"{dataset}: matched controls are not CPM libraries")
        counts = info["nCells"].to_numpy(np.int64)
        implied = truth.to_numpy() * counts[:, None]
        if (
            not np.allclose(truth.sum(axis=1), 1)
            or not np.allclose(implied, np.round(implied), atol=1e-7)
            or not np.array_equal(np.round(implied).sum(axis=1).astype(np.int64), counts)
        ):
            raise ValueError(f"{dataset}: composition denominator failure")

        b = matrix(root / "models/CLAMPfull/B.csv")
        z = matrix(root / "models/CLAMPfull/Z.csv")
        rank = int(pd.read_csv(root / "preprocessing/k.csv")["k"].iloc[0])
        if not np.isfinite(b.to_numpy()).all() or not np.isfinite(z.to_numpy()).all():
            raise ValueError(f"{dataset}: non-finite full model")
        if b.shape != (rank, len(donors)) or z.shape[1] != rank:
            raise ValueError(f"{dataset}: model rank or dimensions are invalid")

        membership = pd.read_csv(root / "grouped_cv/fold_membership.csv", dtype=str)
        membership["fold"] = membership["fold"].astype(int)
        if membership["sample"].duplicated().any() or set(membership["sample"]) != set(donors):
            raise ValueError(f"{dataset}: incomplete fold membership")
        if set(membership.fold) != {1, 2, 3, 4, 5} or membership.groupby("group_id").fold.nunique().max() != 1:
            raise ValueError(f"{dataset}: fold leakage")
        for subdir in ["grouped_cv", "grouped_cv_mean_cell_cpm"]:
            for fold in range(1, 6):
                fold_root = root / subdir / f"fold{fold}/CLAMPfull"
                train = matrix(fold_root / "train_B.csv").T
                test = matrix(fold_root / "test_B.csv").T
                if set(train.index) & set(test.index) or set(train.index) | set(test.index) != set(donors):
                    raise ValueError(f"{dataset}: {subdir} fold {fold} donor coverage failure")
                if not np.isfinite(train.to_numpy()).all() or not np.isfinite(test.to_numpy()).all():
                    raise ValueError(f"{dataset}: {subdir} fold {fold} non-finite scores")

        n_projected, projection_rank = finite_projection(
            root / "single_cell_projection/single_cell_lv_scores.h5"
        )
        projection_summary = pd.read_csv(root / "single_cell_projection/projection_summary.csv").iloc[0]
        if int(projection_summary.n_cells_projected) != n_projected:
            raise ValueError(f"{dataset}: projection coverage failure")
        umap = pd.read_csv(root / "single_cell_umap/umap_points.csv")
        umap_summary = pd.read_csv(root / "single_cell_umap/umap_summary.csv").iloc[0]
        if (
            umap.cell_index.duplicated().any()
            or len(umap) != int(aggregation.n_cells_sampled)
            or len(umap) != min(args.sample_target, int(aggregation.n_cells_retained))
            or not np.isfinite(umap[["umap1", "umap2"]].to_numpy()).all()
            or umap_summary.embedding_source != "sampled_raw_expression_normalize_hvg_pca_neighbors_umap"
        ):
            raise ValueError(f"{dataset}: expression UMAP validation failure")
        rows.append(
            {
                "dataset": dataset,
                "n_genes": donor.shape[0],
                "n_donors": donor.shape[1],
                "n_retained_cells": int(counts.sum()),
                "n_sampled_umap_cells": len(umap),
                "full_rank": rank,
                "n_projected_cells": n_projected,
                "projection_lvs": projection_rank,
                "truth_denominator_exact": True,
                "folds_leakage_free": True,
                "all_model_matrices_finite": True,
                "expression_umap_independent": True,
            }
        )

    metrics = pd.read_csv(analysis / "oof_metrics.csv")
    valid = metrics[metrics.valid_all_folds]
    recovery = pd.read_csv(analysis / "single_cell_recovery.csv")
    observed = {
        "n_types": int(len(metrics)),
        "n_valid": int(len(valid)),
        "mean_oof_r": float(valid.pearson_r.mean()),
        "mean_predictive_r2": float(valid.predictive_r2.mean()),
        "pooled_top1_purity_pct": float(recovery.diag_count.sum() / recovery.n_top.sum() * 100),
    }
    expected = {
        "n_types": 48,
        "n_valid": 32,
        "mean_oof_r": 0.792,
        "mean_predictive_r2": 0.629,
        "pooled_top1_purity_pct": 68.3,
    }
    if observed["n_types"] != expected["n_types"] or observed["n_valid"] != expected["n_valid"]:
        raise ValueError(f"benchmark universe changed: {observed}")
    for key, tolerance in [("mean_oof_r", 0.01), ("mean_predictive_r2", 0.01), ("pooled_top1_purity_pct", 0.5)]:
        if abs(observed[key] - expected[key]) > tolerance:
            raise ValueError(f"benchmark result {key} changed: {observed[key]} vs {expected[key]}")

    frame = pd.DataFrame(rows)
    summary_path, json_path = resolve(args.output_summary), resolve(args.output_json)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    json_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(summary_path, index=False)
    json_path.write_text(
        json.dumps(
            {
                "status": "PASS",
                "datasets": rows,
                "observed_benchmark": observed,
                "expected_benchmark": expected,
            },
            indent=2,
        )
    )
    print(frame.to_string(index=False))
    print("donor-bulk validation PASS")


if __name__ == "__main__":
    main()
