#!/usr/bin/env python3
"""Analyze donor-bulk CLAMPfull recovery, controls, and projected cells."""

from __future__ import annotations

import argparse
import math
from collections import Counter
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.stats import rankdata


REPO_ROOT = Path(__file__).resolve().parents[2]


def resolve(value: str | Path) -> Path:
    value = Path(value)
    return value if value.is_absolute() else REPO_ROOT / value


def read_matrix(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, index_col=0)
    frame.index = frame.index.astype(str)
    frame.columns = frame.columns.astype(str)
    return frame.astype(float)


def safe_cor(x, y) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 2 or np.var(x[keep]) <= 0 or np.var(y[keep]) <= 0:
        return np.nan
    return float(np.corrcoef(x[keep], y[keep])[0, 1])


def correlation_matrix(scores: pd.DataFrame, truth: pd.DataFrame) -> pd.DataFrame:
    shared = sorted(set(scores.index) & set(truth.index))
    scores = scores.loc[shared]
    truth = truth.loc[shared]
    result = np.empty((scores.shape[1], truth.shape[1]), dtype=float)
    for i, lv in enumerate(scores.columns):
        for j, cell_type in enumerate(truth.columns):
            result[i, j] = safe_cor(scores[lv], truth[cell_type])
    return pd.DataFrame(result, index=scores.columns, columns=truth.columns)


def greedy_assign(correlation: pd.DataFrame) -> list[tuple[str, str, float]]:
    work = correlation.to_numpy(copy=True)
    assignments = []
    while np.isfinite(work).any():
        best = np.nanmax(work)
        i, j = np.argwhere(work == best)[0]
        assignments.append((correlation.index[i], correlation.columns[j], float(work[i, j])))
        work[i, :] = np.nan
        work[:, j] = np.nan
    return assignments


def filter_rare(truth: pd.DataFrame) -> pd.DataFrame:
    return truth.loc[:, truth.mean(axis=0) >= 0.005] if len(truth) <= 100 else truth


def regression_predict(x_train, y_train, x_test):
    design = np.column_stack([np.ones(len(x_train)), x_train])
    coefficient, *_ = np.linalg.lstsq(design, y_train, rcond=None)
    return coefficient[0] + coefficient[1] * x_test, coefficient


def full_data_analysis(production: Path, datasets: list[str], output: Path) -> pd.DataFrame:
    assignment_rows, correlation_rows, summary_rows = [], [], []
    for dataset in datasets:
        root = production / dataset
        scores = read_matrix(root / "models/CLAMPfull/B.csv").T
        truth = filter_rare(read_matrix(root / "bulk/truthFrac_v0.csv"))
        correlation = correlation_matrix(scores, truth)
        assignments = greedy_assign(correlation)
        for lv, cell_type, value in assignments:
            other_lv = correlation.loc[correlation.index != lv, cell_type]
            other_type = correlation.loc[lv, correlation.columns != cell_type]
            assignment_rows.append(
                {
                    "dataset": dataset,
                    "cell_type": cell_type,
                    "LV": lv,
                    "cor": value,
                    "r_next_best_lv": float(other_lv.max()) if len(other_lv) else np.nan,
                    "r_next_best_ct": float(other_type.max()) if len(other_type) else np.nan,
                    "margin_lv": value - float(other_lv.max()) if len(other_lv) else np.nan,
                    "margin_ct": value - float(other_type.max()) if len(other_type) else np.nan,
                }
            )
        for lv in correlation.index:
            for cell_type in correlation.columns:
                correlation_rows.append(
                    {
                        "dataset": dataset,
                        "LV": lv,
                        "cell_type": cell_type,
                        "cor": correlation.loc[lv, cell_type],
                    }
                )
        values = [item[2] for item in assignments]
        summary_rows.append(
            {
                "dataset": dataset,
                "n_donors": len(scores),
                "n_lvs": scores.shape[1],
                "n_cell_types": len(truth.columns),
                "n_assigned": len(values),
                "mean_assigned_cor": np.mean(values),
                "median_assigned_cor": np.median(values),
            }
        )
    assignments = pd.DataFrame(assignment_rows)
    assignments.to_csv(output / "full_lv_assignments.csv", index=False)
    pd.DataFrame(correlation_rows).to_csv(output / "full_lv_celltype_correlations.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(output / "full_recovery_by_dataset.csv", index=False)
    return assignments


def grouped_cv_analysis(
    production: Path,
    datasets: list[str],
    output: Path,
    model_subdir: str,
    prefix: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    predictions, calibrations = [], []
    for dataset in datasets:
        root = production / dataset
        truth = filter_rare(read_matrix(root / "bulk/truthFrac_v0.csv"))
        for fold in range(1, 6):
            fold_root = root / model_subdir / f"fold{fold}/CLAMPfull"
            train_b = read_matrix(fold_root / "train_B.csv").T
            test_b = read_matrix(fold_root / "test_B.csv").T
            train_ids = sorted(set(train_b.index) & set(truth.index))
            test_ids = sorted(set(test_b.index) & set(truth.index))
            if set(train_ids) & set(test_ids):
                raise ValueError(f"{dataset} fold {fold}: donor leakage")
            train_truth = truth.loc[train_ids]
            test_truth = truth.loc[test_ids]
            train_b = train_b.loc[train_ids]
            test_b = test_b.loc[test_ids]
            assignment = {
                cell_type: (lv, value)
                for lv, cell_type, value in correlation_assignments(train_b, train_truth)
            }
            for cell_type in train_truth.columns:
                lv, train_r = assignment.get(cell_type, (None, np.nan))
                status = "ok"
                if lv is None or not np.isfinite(train_r):
                    status = "no_one_to_one_lv"
                elif np.var(train_b[lv]) <= 0 or np.var(train_truth[cell_type]) <= 0:
                    status = "constant_training_signal"
                if status == "ok":
                    predicted, coefficient = regression_predict(
                        train_b[lv].to_numpy(),
                        train_truth[cell_type].to_numpy(),
                        test_b[lv].to_numpy(),
                    )
                    alpha, beta = coefficient
                else:
                    alpha, beta = train_truth[cell_type].mean(), 0.0
                    predicted = np.full(len(test_ids), alpha)
                eligible = status == "ok" and train_r >= 0.5
                calibrations.append(
                    {
                        "dataset": dataset,
                        "fold": fold,
                        "cell_type": cell_type,
                        "selected_lv": lv,
                        "train_lv_cor": train_r,
                        "eligible_train_signal": eligible,
                        "alpha": alpha,
                        "beta": beta,
                        "n_train": len(train_ids),
                        "n_test": len(test_ids),
                        "status": status,
                    }
                )
                for sample, observed, estimate in zip(
                    test_ids, test_truth[cell_type], predicted
                ):
                    predictions.append(
                        {
                            "dataset": dataset,
                            "fold": fold,
                            "sample": sample,
                            "cell_type": cell_type,
                            "selected_lv": lv,
                            "observed": observed,
                            "predicted": estimate,
                            "train_lv_cor": train_r,
                            "eligible_train_signal": eligible,
                            "status": status,
                        }
                    )
    predictions = pd.DataFrame(predictions)
    calibrations = pd.DataFrame(calibrations)
    metrics = []
    for (dataset, cell_type), group in predictions.groupby(["dataset", "cell_type"], sort=False):
        folds_valid = group.groupby("fold")["eligible_train_signal"].all().reindex(range(1, 6), fill_value=False)
        observed = group["observed"].to_numpy()
        predicted = group["predicted"].to_numpy()
        sst = np.sum((observed - observed.mean()) ** 2)
        metrics.append(
            {
                "dataset": dataset,
                "cell_type": cell_type,
                "valid_all_folds": bool(folds_valid.all()),
                "n_oof": len(group),
                "pearson_r": safe_cor(predicted, observed),
                "predictive_r2": 1 - np.sum((observed - predicted) ** 2) / sst if sst > 0 else np.nan,
                "mae": np.mean(np.abs(observed - predicted)),
                "minimum_train_lv_cor": group["train_lv_cor"].min(),
            }
        )
    metrics = pd.DataFrame(metrics)
    dataset_rows = []
    for dataset, group in metrics.groupby("dataset", sort=False):
        selected = group[group.valid_all_folds]
        dataset_rows.append(
            {
                "dataset": dataset,
                "n_valid": len(selected),
                "n_total": len(group),
                "mean_pearson_r": selected.pearson_r.mean(),
                "mean_predictive_r2": selected.predictive_r2.mean(),
                "mean_mae": selected.mae.mean(),
            }
        )
    dataset_summary = pd.DataFrame(dataset_rows)
    valid = metrics[metrics.valid_all_folds]
    overall = pd.DataFrame(
        [
            {
                "n_valid": len(valid),
                "n_total": len(metrics),
                "mean_pearson_r": valid.pearson_r.mean(),
                "mean_predictive_r2": valid.predictive_r2.mean(),
                "mean_mae": valid.mae.mean(),
            }
        ]
    )
    predictions.to_csv(output / f"{prefix}oof_predictions.csv", index=False)
    calibrations.to_csv(output / f"{prefix}fold_calibrations.csv", index=False)
    metrics.to_csv(output / f"{prefix}oof_metrics.csv", index=False)
    dataset_summary.to_csv(output / f"{prefix}cv_recovery_by_dataset.csv", index=False)
    overall.to_csv(output / f"{prefix}cv_recovery_overall.csv", index=False)
    return metrics, dataset_summary


def correlation_assignments(scores: pd.DataFrame, truth: pd.DataFrame):
    return greedy_assign(correlation_matrix(scores, truth))


def decode(values) -> np.ndarray:
    return np.asarray(
        [value.decode() if isinstance(value, bytes) else str(value) for value in values],
        dtype=str,
    )


def auc_from_ranks(ranks: np.ndarray, positive: np.ndarray) -> float:
    n_positive = int(positive.sum())
    n_negative = len(positive) - n_positive
    if n_positive == 0 or n_negative == 0:
        return np.nan
    return (
        float(ranks[positive].sum()) - n_positive * (n_positive + 1) / 2
    ) / (n_positive * n_negative)


def single_cell_analysis(
    production: Path,
    datasets: list[str],
    assignments: pd.DataFrame,
    output: Path,
) -> None:
    recovery_rows, specificity_rows, specificity_summary = [], [], []
    for dataset in datasets:
        selected = assignments[assignments.dataset.eq(dataset)].copy()
        truth = read_matrix(production / dataset / "bulk/truthFrac_v0.csv")
        cell_types = [value for value in truth.columns if value in set(selected.cell_type)]
        selected = selected.set_index("cell_type").loc[cell_types].reset_index()
        projection = production / dataset / "single_cell_projection/single_cell_lv_scores.h5"
        with h5py.File(projection, "r") as h5:
            labels = decode(h5["mapped_cell_type"][:])
            lv_names = decode(h5["lv_names"][:])
            lv_lookup = {name: index for index, name in enumerate(lv_names)}
            n_cells = len(labels)
            n_top = max(1, math.ceil(n_cells * 0.01))
            population = Counter(labels)
            positives = {cell_type: labels == cell_type for cell_type in cell_types}
            dataset_specificity = []
            for assigned_order, record in enumerate(selected.itertuples(index=False), start=1):
                scores = h5["scores"][:, lv_lookup[record.LV]].astype(np.float64)
                top = np.argpartition(scores, len(scores) - n_top)[-n_top:]
                diagonal = int(np.sum(labels[top] == record.cell_type))
                n_true = int(population[record.cell_type])
                purity = diagonal / n_top
                prevalence = n_true / n_cells
                maximum = min(1.0, n_true / n_top)
                recovery_rows.append(
                    {
                        "dataset": dataset,
                        "cell_type": record.cell_type,
                        "assigned_lv": record.LV,
                        "assignment_cor": record.cor,
                        "n_top": n_top,
                        "diag_count": diagonal,
                        "recovery_pct": purity * 100,
                        "n_true_cells": n_true,
                        "n_population": n_cells,
                        "prevalence": prevalence,
                        "max_achievable_purity": maximum,
                        "purity_ratio": purity / maximum if maximum else np.nan,
                        "purity_lift": purity / prevalence if prevalence else np.nan,
                    }
                )
                ranks = rankdata(scores, method="average")
                score_mean, score_sd = float(scores.mean()), float(scores.std())
                for observed_order, observed_type in enumerate(cell_types, start=1):
                    positive = positives[observed_type]
                    row = {
                        "dataset": dataset,
                        "assigned_cell_type": record.cell_type,
                        "assigned_lv": record.LV,
                        "observed_cell_type": observed_type,
                        "assigned_order": assigned_order,
                        "observed_order": observed_order,
                        "auc": auc_from_ranks(ranks, positive),
                        "mean_z": (float(scores[positive].mean()) - score_mean) / score_sd
                        if positive.any() and score_sd > 0 else np.nan,
                        "is_match": record.cell_type == observed_type,
                        "n_positive_cells": int(positive.sum()),
                        "n_total_cells": n_cells,
                    }
                    dataset_specificity.append(row)
                    specificity_rows.append(row)
        frame = pd.DataFrame(dataset_specificity)
        diagonal = frame[frame.is_match]
        best = frame.loc[frame.groupby("assigned_cell_type")["auc"].idxmax()]
        specificity_summary.append(
            {
                "dataset": dataset,
                "n_types": len(cell_types),
                "median_matched_auc": diagonal.auc.median(),
                "mean_matched_auc": diagonal.auc.mean(),
                "n_matched_auc_ge_0_7": int((diagonal.auc >= 0.7).sum()),
                "n_matched_is_best": int((best.assigned_cell_type == best.observed_cell_type).sum()),
            }
        )
    recovery = pd.DataFrame(recovery_rows)
    recovery.to_csv(output / "single_cell_recovery.csv", index=False)
    recovery.groupby("dataset", sort=False).agg(
        n_cell_types=("cell_type", "size"),
        mean_recovery_pct=("recovery_pct", "mean"),
        median_recovery_pct=("recovery_pct", "median"),
        mean_purity_ratio=("purity_ratio", "mean"),
        mean_purity_lift=("purity_lift", "mean"),
        median_assignment_cor=("assignment_cor", "median"),
    ).reset_index().to_csv(output / "single_cell_recovery_by_dataset.csv", index=False)
    pd.DataFrame(
        [
            {
                "global_recovery_pct_pooled": recovery.diag_count.sum() / recovery.n_top.sum() * 100,
                "mean_recovery_pct_by_dataset": recovery.groupby("dataset").recovery_pct.mean().mean(),
                "mean_purity_ratio": recovery.purity_ratio.mean(),
                "mean_purity_lift": recovery.purity_lift.mean(),
            }
        ]
    ).to_csv(output / "single_cell_recovery_overall.csv", index=False)
    pd.DataFrame(specificity_rows).to_csv(output / "single_cell_specificity_matrix.csv", index=False)
    pd.DataFrame(specificity_summary).to_csv(output / "single_cell_specificity_summary.csv", index=False)


def comparison(donor: pd.DataFrame, control: pd.DataFrame, output: Path) -> None:
    merged = donor.merge(
        control,
        on=["dataset", "cell_type"],
        suffixes=("_donor_sum", "_mean_cell_cpm"),
    )
    for metric in ["pearson_r", "predictive_r2", "mae"]:
        merged[f"delta_{metric}"] = merged[f"{metric}_donor_sum"] - merged[f"{metric}_mean_cell_cpm"]
    merged.to_csv(output / "comparison_cv_samecell.csv", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--production-root", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--datasets", nargs="+", required=True)
    args = parser.parse_args()
    production = resolve(args.production_root)
    output = resolve(args.output_dir)
    output.mkdir(parents=True, exist_ok=True)
    assignments = full_data_analysis(production, args.datasets, output)
    donor_metrics, _ = grouped_cv_analysis(
        production, args.datasets, output, "grouped_cv", ""
    )
    control_metrics, _ = grouped_cv_analysis(
        production, args.datasets, output, "grouped_cv_mean_cell_cpm", "control_"
    )
    comparison(donor_metrics, control_metrics, output)
    single_cell_analysis(production, args.datasets, assignments, output)
    print("donor-bulk analysis complete")


if __name__ == "__main__":
    main()
