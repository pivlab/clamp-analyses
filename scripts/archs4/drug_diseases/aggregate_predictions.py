#!/usr/bin/env python3
# Aggregate per-tissue drug-disease prediction HDF5 files into report tables.
# Loads every method's HDF5 prediction files, ranks scores within the full DOID
# distribution, merges with the gold standard, reduces across thresholds (mean)
# then across tissues (max), and writes performance_summary.csv (AUROC/AvgPrecision
# per method) plus the raw/aggregated prediction tables consumed by notebooks 00/01.

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

N_TISSUES = 49

METHOD_RENAME = {
    "Gene-based": "gene_based",
    "Module-based (ARCHS4)": "module_based_archs4",
    "Module-based": "module_based_archs4",
    "Module-based (GTEx)": "module_based_gtex",
    "Module-based (recount2)": "module_based_recount2",
}

def parse_thresholds(raw: str) -> list[float]:
    return [-1.0 if tok == "all" else float(tok) for tok in raw.split(",")]


def _get_tissue(data_value: str) -> str:
    prefix = "spredixcan-mashr-zscores-"
    assert data_value.startswith(prefix), data_value
    tissue = data_value[len(prefix):]
    for suffix in (
        "-projection-archs4",
        "-projection-gtex",
        "-projection-recount2",
        "-projection",
        "-data",
    ):
        if tissue.endswith(suffix):
            return tissue[: -len(suffix)]
    raise ValueError(f"Cannot extract tissue from metadata data value: {data_value}")


def load_predictions(
    prediction_dirs: dict[str, Path],
    gold_standard: pd.DataFrame,
    method_thresholds: dict[str, list[float]],
) -> pd.DataFrame:
    expected_methods = tuple(method_thresholds)

    current_prediction_files = []
    for d in prediction_dirs.values():
        current_prediction_files.extend(sorted(d.glob("*.h5")))
    current_prediction_files.sort()
    print(f"Found {len(current_prediction_files)} prediction files")

    predictions = []
    skipped_files = []
    for f in current_prediction_files:
        metadata = pd.read_hdf(f, key="metadata")
        method_name = METHOD_RENAME.get(metadata["method"].values[0], metadata["method"].values[0])
        if method_name not in method_thresholds:
            skipped_files.append((f.name, method_name))
            continue

        prediction_data = pd.read_hdf(f, key="prediction")
        prediction_data["score"] = prediction_data["score"].rank()
        prediction_data = pd.merge(prediction_data, gold_standard, on=["trait", "drug"], how="inner")
        prediction_data["trait"] = prediction_data["trait"].astype("category")
        prediction_data["drug"] = prediction_data["drug"].astype("category")

        prediction_data = prediction_data.assign(method=method_name)
        prediction_data["method"] = pd.Categorical(
            prediction_data["method"], categories=expected_methods, ordered=True
        )
        prediction_data = prediction_data.assign(n_top_genes=metadata["n_top_genes"].values[0])

        data_value = metadata["data"].values[0]
        prediction_data = prediction_data.assign(data=data_value)
        prediction_data["data"] = prediction_data["data"].astype("category")
        prediction_data = prediction_data.assign(tissue=_get_tissue(data_value))

        predictions.append(prediction_data)

    print(f"Skipped files: {len(skipped_files)}")
    if skipped_files:
        print(skipped_files[:10])

    predictions = pd.concat(predictions, ignore_index=True)
    assert not predictions.isna().any().any()

    method_counts = predictions["method"].value_counts().reindex(expected_methods)
    n_predictions = predictions[["drug", "trait"]].drop_duplicates().shape[0]
    print(f"Unique drug-disease pairs: {n_predictions}")

    for method_name, thresholds in method_thresholds.items():
        expected = N_TISSUES * len(thresholds) * n_predictions
        actual = int(method_counts.loc[method_name])
        assert actual == expected, f"{method_name}: expected {expected}, got {actual}"

    threshold_counts = predictions.groupby(["method", "n_top_genes"], observed=True).size().rename("n")
    for method_name, thresholds in method_thresholds.items():
        counts = threshold_counts.loc[method_name]
        actual_thresholds = sorted(float(x) for x in counts.index)
        expected_thresholds = sorted(thresholds)
        assert actual_thresholds == expected_thresholds, (
            f"{method_name}: expected thresholds {expected_thresholds}, got {actual_thresholds}"
        )
        assert np.all(counts.values == N_TISSUES * n_predictions), method_name

    return predictions


def _reduce_mean(x: pd.DataFrame) -> pd.Series:
    return pd.Series({"score": x["score"].mean(), "true_class": x["true_class"].unique()[0]})


def _reduce_max(x: pd.DataFrame) -> pd.Series:
    return pd.Series({"score": x["score"].max(), "true_class": x["true_class"].unique()[0]})


def reduce_mean_then_max(predictions: pd.DataFrame) -> pd.DataFrame:
    return (
        predictions
        .groupby(["trait", "drug", "method", "tissue"], observed=True)
        .apply(_reduce_mean, include_groups=False)
        .dropna()
        .groupby(["trait", "drug", "method"], observed=True)
        .apply(_reduce_max, include_groups=False)
        .dropna()
        .sort_index()
        .reset_index()
    )


def build_method_thresholds(n_top_lvs: list[float], n_top_genes: list[float]) -> dict[str, list[float]]:
    return {
        "module_based_archs4": n_top_lvs,
        "module_based_gtex": n_top_lvs,
        "module_based_recount2": n_top_lvs,
        "gene_based": n_top_genes,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--module-archs4-dir", required=True, type=Path)
    parser.add_argument("--module-gtex-dir", required=True, type=Path)
    parser.add_argument("--module-recount2-dir", required=True, type=Path)
    parser.add_argument("--gene-dir", required=True, type=Path)
    parser.add_argument("--gold-standard", required=True, type=Path)
    parser.add_argument("--ukb-efo", required=True, type=Path)
    parser.add_argument("--efo-xrefs", required=True, type=Path)
    parser.add_argument("--do-xrefs", required=True, type=Path)
    parser.add_argument("--n-top-lvs", required=True, help="comma-separated ints, or 'all' entries")
    parser.add_argument("--n-top-genes", required=True, help="comma-separated ints, or 'all' entries")
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    prediction_dirs = {
        "module_based_archs4": args.module_archs4_dir / "lincs" / "predictions" / "dotprod_neg",
        "module_based_gtex": args.module_gtex_dir / "lincs" / "predictions" / "dotprod_neg",
        "module_based_recount2": args.module_recount2_dir / "lincs" / "predictions" / "dotprod_neg",
        "gene_based": args.gene_dir / "lincs" / "predictions" / "dotprod_neg",
    }
    for name, d in prediction_dirs.items():
        assert d.exists(), f"{name} predictions missing: {d}"

    method_thresholds = build_method_thresholds(parse_thresholds(args.n_top_lvs), parse_thresholds(args.n_top_genes))

    gold_standard = pd.read_pickle(args.gold_standard)
    print(f"Gold standard shape: {gold_standard.shape}")

    predictions = load_predictions(prediction_dirs, gold_standard, method_thresholds)
    print(f"Loaded predictions shape: {predictions.shape}")

    predictions.to_pickle(args.output_dir / "predictions_results.pkl")

    predictions_avg = reduce_mean_then_max(predictions)
    assert predictions_avg.shape[0] == len(method_thresholds) * predictions[["drug", "trait"]].drop_duplicates().shape[0]
    predictions_avg.to_pickle(args.output_dir / "predictions_results_aggregated.pkl")

    auroc_final = predictions_avg.groupby("method", observed=True).apply(
        lambda x: roc_auc_score(x["true_class"], x["score"]), include_groups=False
    ).rename("AUROC")
    ap_final = predictions_avg.groupby("method", observed=True).apply(
        lambda x: average_precision_score(x["true_class"], x["score"]), include_groups=False
    ).rename("AvgPrecision")
    summary = pd.concat([auroc_final, ap_final], axis=1)
    print(summary)
    summary.to_csv(args.output_dir / "performance_summary.csv")


if __name__ == "__main__":
    main()
