#!/usr/bin/env python3
# Evaluate donor-grouped GTEx anatomical subtissue recovery.
# test recovery of detailed SMTSD labels within each broad tissue.


from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import rpy2.robjects as ro
from joblib import Parallel, delayed, parallel_config
from sklearn.dummy import DummyClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import extract_B_matrix  # noqa: E402

readRDS = ro.r["readRDS"]

REPRESENTATIONS = (
    ("FullLV", "Full CLAMP LV space"),
    ("SHAP", "Broad-tissue RF SHAP"),
)
META_COLUMNS = {"SAMPID", "donor_id", "fold", "Tissue", "SMTSD"}


def _short_label(label: str, tissue: str) -> str:
    prefix = f"{tissue} - "
    if label.startswith(prefix):
        return label[len(prefix):]
    return label


def _fit_predict_lr(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: np.ndarray,
    seed: int,
) -> np.ndarray:
    if len(np.unique(y_train)) < 2:
        return np.full(len(X_test), y_train[0], dtype=object)
    classifier = LogisticRegression(
        max_iter=2_000,
        class_weight="balanced",
        solver="lbfgs",
        random_state=seed,
    )
    classifier.fit(X_train, y_train)
    return classifier.predict(X_test)


def _fit_predict_chance(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: np.ndarray,
    seed: int,
) -> np.ndarray:
    classifier = DummyClassifier(strategy="stratified", random_state=seed)
    classifier.fit(X_train, y_train)
    return classifier.predict(X_test)


def _oof_predict(
    features: np.ndarray,
    labels: np.ndarray,
    folds: np.ndarray,
    seed: int,
    predict_fn=_fit_predict_lr,
) -> np.ndarray:
    predictions = np.empty(len(labels), dtype=object)
    for fold in np.sort(np.unique(folds)):
        train_mask = folds != fold
        test_mask = folds == fold
        predictions[test_mask] = predict_fn(
            features[train_mask], labels[train_mask], features[test_mask], seed + int(fold),
        )
    assert not pd.isna(predictions).any(), "OOF prediction vector contains missing values"
    return predictions


def _balanced_accuracy_fixed(
    labels: np.ndarray,
    predictions: np.ndarray,
    classes: np.ndarray,
) -> float:
    recalls = []
    for label in classes:
        mask = labels == label
        if not mask.any():
            raise ValueError(f"class missing while calculating balanced accuracy: {label}")
        recalls.append(float(np.mean(predictions[mask] == label)))
    return float(np.mean(recalls))


def _metrics(
    labels: np.ndarray,
    predictions: np.ndarray,
    classes: np.ndarray,
) -> dict[str, float]:
    balanced_accuracy = _balanced_accuracy_fixed(labels, predictions, classes)
    sklearn_balanced_accuracy = balanced_accuracy_score(labels, predictions)
    assert np.isclose(balanced_accuracy, sklearn_balanced_accuracy), (
        "balanced accuracy does not equal mean per-class recall"
    )
    n_classes = len(classes)
    adjusted = (balanced_accuracy - 1.0 / n_classes) / (1.0 - 1.0 / n_classes)
    return {
        "Balanced_Accuracy": balanced_accuracy,
        "Adjusted_Balanced_Accuracy": adjusted,
        "Macro_F1": float(
            f1_score(labels, predictions, labels=classes, average="macro", zero_division=0)
        ),
        "Accuracy": float(accuracy_score(labels, predictions)),
    }


def _per_class_recall(
    labels: np.ndarray,
    predictions: np.ndarray,
    classes: np.ndarray,
) -> dict[str, float]:
    return {
        label: float(np.mean(predictions[labels == label] == label))
        for label in classes
    }


def _permute_donor_profiles(
    labels: np.ndarray,
    donors: np.ndarray,
    rng: np.random.Generator,
) -> np.ndarray:
    """Permute complete donor label profiles within donor-size strata.

    Donor blocks remain intact.  Assigning each source profile exactly once
    preserves global class counts, while random within-profile alignment
    breaks the sample-level association between representation and SMTSD.
    """
    donor_indices = {
        donor: np.flatnonzero(donors == donor)
        for donor in np.unique(donors)
    }
    donors_by_size: dict[int, list[str]] = defaultdict(list)
    for donor, indices in donor_indices.items():
        donors_by_size[len(indices)].append(donor)

    permuted = np.empty_like(labels, dtype=object)
    for donor_size, target_donors in sorted(donors_by_size.items()):
        target_donors = sorted(target_donors)
        source_donors = np.asarray(target_donors, dtype=object)[rng.permutation(len(target_donors))]
        for target_donor, source_donor in zip(target_donors, source_donors):
            target_indices = donor_indices[target_donor]
            source_labels = labels[donor_indices[source_donor]]
            assert len(source_labels) == donor_size
            permuted[target_indices] = rng.permutation(source_labels)

    assert Counter(permuted) == Counter(labels), "donor-aware permutation changed class counts"
    return permuted


def _permutation_worker(
    permutation: int,
    features_by_representation: dict[str, np.ndarray],
    labels: np.ndarray,
    donors: np.ndarray,
    folds: np.ndarray,
    classes: np.ndarray,
    seed: int,
) -> list[dict[str, object]]:
    rng = np.random.default_rng(np.random.SeedSequence([seed, permutation]))
    permuted_labels = _permute_donor_profiles(labels, donors, rng)
    rows = []
    for representation, features in features_by_representation.items():
        predictions = _oof_predict(
            features, permuted_labels, folds, seed + permutation * 10, _fit_predict_lr,
        )
        rows.append({
            "Permutation": permutation,
            "Representation": representation,
            "Balanced_Accuracy": _balanced_accuracy_fixed(
                permuted_labels, predictions, classes,
            ),
        })
    return rows


def _bootstrap_predictions(
    labels: np.ndarray,
    donors: np.ndarray,
    predictions_by_representation: dict[str, np.ndarray],
    classes: np.ndarray,
    n_repeats: int,
    seed: int,
) -> list[dict[str, object]]:
    rng = np.random.default_rng(seed)
    unique_donors = np.sort(np.unique(donors))
    donor_indices = {
        donor: np.flatnonzero(donors == donor)
        for donor in unique_donors
    }
    rows = []
    for bootstrap in range(1, n_repeats + 1):
        for _attempt in range(100):
            sampled_donors = rng.choice(unique_donors, size=len(unique_donors), replace=True)
            sampled_indices = np.concatenate([donor_indices[donor] for donor in sampled_donors])
            sampled_labels = labels[sampled_indices]
            if set(np.unique(sampled_labels)) == set(classes):
                break
        else:
            raise RuntimeError("unable to draw a donor bootstrap containing every class")

        for representation, predictions in predictions_by_representation.items():
            score = _balanced_accuracy_fixed(
                sampled_labels, predictions[sampled_indices], classes,
            )
            rows.append({
                "Bootstrap": bootstrap,
                "Representation": representation,
                "Balanced_Accuracy": score,
                "Adjusted_Balanced_Accuracy": (
                    (score - 1.0 / len(classes)) / (1.0 - 1.0 / len(classes))
                ),
            })
    return rows


def _complete_confusion_rows(
    labels: np.ndarray,
    predictions: np.ndarray,
    classes: np.ndarray,
    class_display: dict[str, str],
) -> list[dict[str, object]]:
    rows = []
    for true_label in classes:
        true_mask = labels == true_label
        denominator = int(true_mask.sum())
        for predicted_label in classes:
            count = int(np.sum(true_mask & (predictions == predicted_label)))
            rows.append({
                "True_SMTSD": true_label,
                "True_Display": class_display[true_label],
                "Predicted_SMTSD": predicted_label,
                "Predicted_Display": class_display[predicted_label],
                "Count": count,
                "Row_Proportion": count / denominator,
            })
    row_sums = pd.DataFrame(rows).groupby("True_SMTSD")["Row_Proportion"].sum()
    assert np.allclose(row_sums, 1.0), "row-normalized confusion matrix does not sum to one"
    return rows


def _evaluate_task(
    tissue: str,
    task_labels: list[str],
    tissue_shap: pd.DataFrame,
    lv_matrix: pd.DataFrame,
    lv_columns: list[str],
    class_display: dict[str, str],
    n_folds: int,
    n_permutations: int,
    n_bootstraps: int,
    n_jobs: int,
    seed: int,
) -> dict[str, object]:
    data = tissue_shap[tissue_shap["SMTSD"].isin(task_labels)].copy()
    data = data.sort_values("SAMPID").reset_index(drop=True)
    classes = np.asarray(task_labels, dtype=object)
    assert set(data["SMTSD"]) == set(classes)
    assert data["SAMPID"].is_unique
    assert set(data["fold"].unique()) == set(range(1, n_folds + 1))
    per_class_folds = data.groupby("SMTSD")["fold"].nunique()
    assert (per_class_folds == n_folds).all(), f"{tissue}: a class is absent from at least one fold"

    sample_ids = data["SAMPID"].to_numpy()
    labels = data["SMTSD"].to_numpy(dtype=object)
    donors = data["donor_id"].to_numpy(dtype=object)
    folds = data["fold"].to_numpy(dtype=int)
    features_by_representation = {
        "FullLV": lv_matrix.loc[sample_ids, lv_columns].to_numpy(dtype=float),
        "SHAP": data[lv_columns].to_numpy(dtype=float),
    }
    assert features_by_representation["FullLV"].shape == features_by_representation["SHAP"].shape

    predictions_by_representation = {
        representation: _oof_predict(features, labels, folds, seed, _fit_predict_lr)
        for representation, features in features_by_representation.items()
    }
    chance_predictions = _oof_predict(
        features_by_representation["SHAP"], labels, folds, seed, _fit_predict_chance,
    )

    observed_metrics = {
        representation: _metrics(labels, predictions, classes)
        for representation, predictions in predictions_by_representation.items()
    }
    recalls = {
        representation: _per_class_recall(labels, predictions, classes)
        for representation, predictions in predictions_by_representation.items()
    }

    bootstrap_rows = _bootstrap_predictions(
        labels,
        donors,
        predictions_by_representation,
        classes,
        n_bootstraps,
        seed + 20_000,
    )

    print(
        f"[gtex] {tissue} ({len(classes)} classes): running {n_permutations} paired "
        f"donor-aware permutations with n_jobs={n_jobs}",
        flush=True,
    )
    with parallel_config(backend="loky", inner_max_num_threads=1):
        permutation_batches = Parallel(n_jobs=n_jobs, max_nbytes="20M", verbose=5)(
            delayed(_permutation_worker)(
                permutation,
                features_by_representation,
                labels,
                donors,
                folds,
                classes,
                seed,
            )
            for permutation in range(1, n_permutations + 1)
        )
    permutation_rows = [row for batch in permutation_batches for row in batch]

    for representation, _display in REPRESENTATIONS:
        null_scores = np.asarray([
            row["Balanced_Accuracy"]
            for row in permutation_rows
            if row["Representation"] == representation
        ])
        assert len(null_scores) == n_permutations
        observed = observed_metrics[representation]["Balanced_Accuracy"]
        exceedances = int(np.sum(null_scores >= observed))
        p_value = (exceedances + 1) / (n_permutations + 1)
        observed_metrics[representation].update({
            "Permutation_P_Value": p_value,
            "Permutation_Exceedances": exceedances,
            "N_Permutations": n_permutations,
            "P_Value_At_Floor": exceedances == 0,
        })

        bootstrap_scores = np.asarray([
            row["Balanced_Accuracy"]
            for row in bootstrap_rows
            if row["Representation"] == representation
        ])
        ci_lower, ci_upper = np.quantile(bootstrap_scores, [0.025, 0.975])
        n_classes = len(classes)
        observed_metrics[representation].update({
            "BA_CI_Lower": float(ci_lower),
            "BA_CI_Upper": float(ci_upper),
            "Adjusted_BA_CI_Lower": float(
                (ci_lower - 1.0 / n_classes) / (1.0 - 1.0 / n_classes)
            ),
            "Adjusted_BA_CI_Upper": float(
                (ci_upper - 1.0 / n_classes) / (1.0 - 1.0 / n_classes)
            ),
        })

    print(
        f"[gtex] {tissue}: FullLV BA={observed_metrics['FullLV']['Balanced_Accuracy']:.4f}; "
        f"SHAP BA={observed_metrics['SHAP']['Balanced_Accuracy']:.4f}",
        flush=True,
    )
    return {
        "Tissue": tissue,
        "Task_Labels": task_labels,
        "Classes": classes,
        "Sample_IDs": sample_ids,
        "Labels": labels,
        "Donors": donors,
        "Folds": folds,
        "Predictions": predictions_by_representation,
        "Chance_Predictions": chance_predictions,
        "Metrics": observed_metrics,
        "Recalls": recalls,
        "Permutation_Rows": permutation_rows,
        "Bootstrap_Rows": bootstrap_rows,
        "Class_Display": class_display,
        "N_Samples": len(data),
        "N_Donors": int(data["donor_id"].nunique()),
    }


def _wide_result_row(
    result: dict[str, object],
    tissue_display: str,
    tissue_order: int,
    rf_accuracy: pd.DataFrame,
) -> dict[str, object]:
    tissue = result["Tissue"]
    tissue_rf = rf_accuracy[rf_accuracy["Tissue"] == tissue]
    row = {
        "Tissue": tissue,
        "Tissue_Display": tissue_display,
        "Tissue_Order": tissue_order,
        "N_Subtissue_Classes": len(result["Classes"]),
        "N_Samples": result["N_Samples"],
        "N_Donors": result["N_Donors"],
        "RF_Mean_OOF_Balanced_Accuracy": float(tissue_rf["OOF_Balanced_Accuracy"].mean()),
        "RF_Min_OOF_Balanced_Accuracy": float(tissue_rf["OOF_Balanced_Accuracy"].min()),
        "RF_Passed_Accuracy_Gate": bool(tissue_rf["Passed_Min_Accuracy"].all()),
    }
    for representation, _display in REPRESENTATIONS:
        for metric, value in result["Metrics"][representation].items():
            row[f"{representation}_{metric}"] = value
    return row


def _decorate_long_rows(
    rows: list[dict[str, object]],
    analysis: str,
    tissue: str,
    tissue_display: str,
) -> list[dict[str, object]]:
    return [
        {
            "Analysis": analysis,
            "Tissue": tissue,
            "Tissue_Display": tissue_display,
            **row,
        }
        for row in rows
    ]


def _write_legacy_outputs(
    output_dir: Path,
    anatomical_results: pd.DataFrame,
    all_results_by_tissue: dict[str, dict[str, object]],
    rf_accuracy: pd.DataFrame,
    dropped_rows: list[dict[str, object]],
) -> None:
    legacy_summary_rows = []
    for row in anatomical_results.to_dict("records"):
        chance = 1.0 / row["N_Subtissue_Classes"]
        legacy_summary_rows.append({
            "Tissue": row["Tissue"],
            "N_Subtissue_Classes": row["N_Subtissue_Classes"],
            "RF_Mean_OOF_Balanced_Accuracy": row["RF_Mean_OOF_Balanced_Accuracy"],
            "RF_Min_OOF_Balanced_Accuracy": row["RF_Min_OOF_Balanced_Accuracy"],
            "RF_Passed_Accuracy_Gate": row["RF_Passed_Accuracy_Gate"],
            "SHAP_LR_Balanced_Accuracy": row["SHAP_Balanced_Accuracy"],
            "SHAP_LR_Accuracy": row["SHAP_Accuracy"],
            "SHAP_LR_F1_Macro": row["SHAP_Macro_F1"],
            "FullLV_LR_Balanced_Accuracy": row["FullLV_Balanced_Accuracy"],
            "FullLV_LR_Accuracy": row["FullLV_Accuracy"],
            "FullLV_LR_F1_Macro": row["FullLV_Macro_F1"],
            "Chance_Balanced_Accuracy": chance,
            "Chance_Accuracy": chance,
            "Chance_F1_Macro": chance,
            "SHAP_vs_Chance_Diff": row["SHAP_Balanced_Accuracy"] - chance,
            "FullLV_vs_Chance_Diff": row["FullLV_Balanced_Accuracy"] - chance,
            "SHAP_vs_FullLV_Diff": (
                row["SHAP_Balanced_Accuracy"] - row["FullLV_Balanced_Accuracy"]
            ),
            "Permutation_P_Value": row["SHAP_Permutation_P_Value"],
        })
    pd.DataFrame(legacy_summary_rows).to_csv(
        output_dir / "subtissue_lr_summary.tsv", sep="\t", index=False,
    )

    legacy_prediction_rows = []
    legacy_confusion_rows = []
    for tissue, result in all_results_by_tissue.items():
        for index, sample_id in enumerate(result["Sample_IDs"]):
            legacy_prediction_rows.append({
                "Tissue": tissue,
                "SAMPID": sample_id,
                "donor_id": result["Donors"][index],
                "Fold": result["Folds"][index],
                "True_SMTSD": result["Labels"][index],
                "SHAP_LR_Pred": result["Predictions"]["SHAP"][index],
                "FullLV_LR_Pred": result["Predictions"]["FullLV"][index],
                "Chance_Pred": result["Chance_Predictions"][index],
            })
        for representation, legacy_name in (("SHAP", "SHAP_LR"), ("FullLV", "FullLV_LR")):
            for row in _complete_confusion_rows(
                result["Labels"],
                result["Predictions"][representation],
                result["Classes"],
                result["Class_Display"],
            ):
                legacy_confusion_rows.append({
                    "Tissue": tissue,
                    "Model": legacy_name,
                    "True_SMTSD": row["True_SMTSD"],
                    "Predicted_SMTSD": row["Predicted_SMTSD"],
                    "Count": row["Count"],
                })
    pd.DataFrame(legacy_prediction_rows).to_csv(
        output_dir / "subtissue_lr_oof_predictions.tsv", sep="\t", index=False,
    )
    pd.DataFrame(legacy_confusion_rows).to_csv(
        output_dir / "subtissue_lr_confusion_matrices.tsv", sep="\t", index=False,
    )
    pd.DataFrame(dropped_rows, columns=["Tissue", "SMTSD", "N_Samples", "Reason"]).to_csv(
        output_dir / "dropped_subtissues.tsv", sep="\t", index=False,
    )
    rf_accuracy.to_csv(output_dir / "rf_oof_accuracy_summary.tsv", sep="\t", index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--membership", required=True)
    parser.add_argument("--fold-dir", required=True, help="Directory containing fold1..foldN")
    parser.add_argument("--n-folds", type=int, default=5)
    parser.add_argument("--clamp-full-rds", required=True)
    parser.add_argument("--gtex-meta", required=True)
    parser.add_argument("--anatomical-spec", required=True)
    parser.add_argument("--tissues", nargs="+", required=True)
    parser.add_argument("--min-subtissue-samples", type=int, default=20)
    parser.add_argument("--min-rf-accuracy", type=float, default=0.85)
    parser.add_argument("--global-seed", type=int, default=42)
    parser.add_argument("--permutation-repeats", type=int, default=1_000)
    parser.add_argument("--bootstrap-repeats", type=int, default=1_000)
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    output_dir = Path(args.out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    membership = pd.read_csv(args.membership, sep="\t", dtype={"fold": int})
    assert membership["SAMPID"].is_unique
    assert set(membership["fold"].unique()) == set(range(1, args.n_folds + 1))
    assert (membership.groupby("donor_id")["fold"].nunique() == 1).all(), (
        "a donor spans multiple folds"
    )

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
    anatomical_spec = pd.read_csv(args.anatomical_spec, sep="\t")
    required_spec_columns = {
        "Tissue", "SMTSD", "Tissue_Display", "Class_Display", "Tissue_Order", "Class_Order",
    }
    assert required_spec_columns.issubset(anatomical_spec.columns)
    assert not anatomical_spec.duplicated(["Tissue", "SMTSD"]).any()
    assert set(anatomical_spec["SMTSD"]).issubset(set(gtex_meta["SMTSD"])), (
        "anatomical specification contains labels absent from GTEx metadata"
    )
    assert set(anatomical_spec["SMTSD"]).issubset(set(membership["SMTSD"])), (
        "anatomical specification contains labels absent from fold membership"
    )

    expected_anatomical_tissues = {
        "Adipose Tissue", "Blood Vessel", "Brain", "Colon", "Esophagus", "Heart", "Skin",
    }
    assert set(anatomical_spec["Tissue"]) == expected_anatomical_tissues
    assert len(anatomical_spec[anatomical_spec["Tissue"] == "Brain"]) == 13
    assert "Blood" not in set(anatomical_spec["Tissue"])
    assert "Kidney" not in set(anatomical_spec["Tissue"])
    assert "Cells - Cultured fibroblasts" not in set(anatomical_spec["SMTSD"])

    clamp_full = readRDS(str(args.clamp_full_rds))
    lv_matrix = extract_B_matrix(clamp_full).T
    lv_matrix = lv_matrix.loc[membership["SAMPID"]]

    fold_dir = Path(args.fold_dir)
    oof_shap_frames = []
    accuracy_frames = []
    for fold in range(1, args.n_folds + 1):
        oof_shap_frames.append(pd.read_csv(fold_dir / f"fold{fold}" / "oof_shap.tsv", sep="\t"))
        accuracy_frames.append(
            pd.read_csv(fold_dir / f"fold{fold}" / "accuracy_summary.tsv", sep="\t")
        )
    oof_shap = pd.concat(oof_shap_frames, ignore_index=True)
    rf_accuracy = pd.concat(accuracy_frames, ignore_index=True)
    assert not oof_shap.duplicated(["SAMPID", "Tissue"]).any()

    membership_lookup = membership.set_index("SAMPID")
    aligned_membership = membership_lookup.loc[oof_shap["SAMPID"]]
    assert np.array_equal(aligned_membership["donor_id"].to_numpy(), oof_shap["donor_id"].to_numpy())
    assert np.array_equal(aligned_membership["fold"].to_numpy(), oof_shap["fold"].to_numpy())
    assert np.array_equal(aligned_membership["SMTSD"].to_numpy(), oof_shap["SMTSD"].to_numpy())

    lv_columns = [column for column in oof_shap.columns if column not in META_COLUMNS]
    assert len(lv_columns) == lv_matrix.shape[1], "SHAP and full-LV feature counts differ"
    assert set(lv_columns) == set(lv_matrix.columns), "SHAP and full-LV feature names differ"
    lv_matrix = lv_matrix.loc[:, lv_columns]

    configured_membership = membership[membership["SMTS"].isin(args.tissues)].copy()
    detailed_counts = (
        configured_membership.groupby(["SMTS", "SMTSD"])
        .agg(N_Samples=("SAMPID", "size"), N_Donors=("donor_id", "nunique"))
        .reset_index()
    )
    eligible_counts = detailed_counts[
        detailed_counts["N_Samples"] >= args.min_subtissue_samples
    ]
    all_labels_by_tissue = {
        tissue: sorted(group["SMTSD"].tolist())
        for tissue, group in eligible_counts.groupby("SMTS")
        if len(group) >= 2
    }
    assert set(all_labels_by_tissue) == {
        "Adipose Tissue", "Blood", "Blood Vessel", "Brain", "Colon", "Esophagus", "Heart", "Skin",
    }

    anatomical_labels_by_tissue = {
        tissue: group.sort_values("Class_Order")["SMTSD"].tolist()
        for tissue, group in anatomical_spec.groupby("Tissue", sort=False)
    }
    for tissue, labels in anatomical_labels_by_tissue.items():
        assert set(labels).issubset(set(all_labels_by_tissue[tissue]))

    anatomical_lookup = anatomical_spec.set_index(["Tissue", "SMTSD"])
    class_display_by_tissue = {}
    for tissue in args.tissues:
        labels = detailed_counts.loc[detailed_counts["SMTS"] == tissue, "SMTSD"].tolist()
        class_display_by_tissue[tissue] = {
            label: (
                anatomical_lookup.loc[(tissue, label), "Class_Display"]
                if (tissue, label) in anatomical_lookup.index
                else _short_label(label, tissue)
            )
            for label in labels
        }

    # Fold audit includes explicit zero rows for sparse/excluded labels.
    audit_rows = []
    anatomical_pairs = set(map(tuple, anatomical_spec[["Tissue", "SMTSD"]].to_numpy()))
    all_pairs = {
        (tissue, label)
        for tissue, labels in all_labels_by_tissue.items()
        for label in labels
    }
    for row in detailed_counts.itertuples(index=False):
        for fold in range(1, args.n_folds + 1):
            sub = configured_membership[
                (configured_membership["SMTS"] == row.SMTS)
                & (configured_membership["SMTSD"] == row.SMTSD)
                & (configured_membership["fold"] == fold)
            ]
            audit_rows.append({
                "Fold": fold,
                "Tissue": row.SMTS,
                "SMTSD": row.SMTSD,
                "N_Samples": len(sub),
                "N_Donors": sub["donor_id"].nunique(),
                "Included_Anatomical": (row.SMTS, row.SMTSD) in anatomical_pairs,
                "Included_All_SMTSD": (row.SMTS, row.SMTSD) in all_pairs,
            })
    fold_audit = pd.DataFrame(audit_rows)
    assert (fold_audit[fold_audit["Included_All_SMTSD"]]["N_Samples"] > 0).all()
    fold_audit.to_csv(output_dir / "fold_smtsd_audit.tsv", sep="\t", index=False)

    # Evaluate each unique (tissue, label-set) task once; common anatomical/all
    # tasks are reused byte-for-byte in both result sets.
    task_definitions = []
    seen_task_keys = set()
    for tissue, labels in anatomical_labels_by_tissue.items():
        key = (tissue, frozenset(labels))
        task_definitions.append((key, labels))
        seen_task_keys.add(key)
    for tissue, labels in all_labels_by_tissue.items():
        key = (tissue, frozenset(labels))
        if key not in seen_task_keys:
            task_definitions.append((key, labels))
            seen_task_keys.add(key)

    task_results = {}
    for task_index, (task_key, labels) in enumerate(task_definitions):
        tissue = task_key[0]
        tissue_shap = oof_shap[oof_shap["Tissue"] == tissue].copy()
        task_seed = args.global_seed + task_index * 100_000
        task_results[task_key] = _evaluate_task(
            tissue=tissue,
            task_labels=labels,
            tissue_shap=tissue_shap,
            lv_matrix=lv_matrix,
            lv_columns=lv_columns,
            class_display={label: class_display_by_tissue[tissue][label] for label in labels},
            n_folds=args.n_folds,
            n_permutations=args.permutation_repeats,
            n_bootstraps=args.bootstrap_repeats,
            n_jobs=args.n_jobs,
            seed=task_seed,
        )

    anatomical_rows = []
    all_rows = []
    analysis_task_maps = {"anatomical": {}, "all_smtsd": {}}
    for tissue, labels in anatomical_labels_by_tissue.items():
        result = task_results[(tissue, frozenset(labels))]
        spec_group = anatomical_spec[anatomical_spec["Tissue"] == tissue]
        anatomical_rows.append(_wide_result_row(
            result,
            spec_group["Tissue_Display"].iloc[0],
            int(spec_group["Tissue_Order"].iloc[0]),
            rf_accuracy,
        ))
        analysis_task_maps["anatomical"][tissue] = result

    all_order = {
        tissue: index + 1
        for index, tissue in enumerate([
            "Adipose Tissue", "Blood Vessel", "Brain", "Colon", "Esophagus", "Heart", "Skin", "Blood",
        ])
    }
    for tissue, labels in all_labels_by_tissue.items():
        result = task_results[(tissue, frozenset(labels))]
        all_rows.append(_wide_result_row(
            result,
            f"{tissue} ({len(labels)} detailed labels)",
            all_order[tissue],
            rf_accuracy,
        ))
        analysis_task_maps["all_smtsd"][tissue] = result

    anatomical_results = pd.DataFrame(anatomical_rows).sort_values("Tissue_Order").reset_index(drop=True)
    all_results = pd.DataFrame(all_rows).sort_values("Tissue_Order").reset_index(drop=True)
    assert len(anatomical_results) == 7
    assert len(all_results) == 8
    anatomical_results.to_csv(
        output_dir / "anatomical_subtissue_results.tsv", sep="\t", index=False,
    )
    all_results.to_csv(output_dir / "all_smtsd_results.tsv", sep="\t", index=False)

    prediction_rows = []
    confusion_rows = []
    permutation_rows = []
    bootstrap_rows = []
    for analysis, tissue_results in analysis_task_maps.items():
        result_table = anatomical_results if analysis == "anatomical" else all_results
        display_lookup = result_table.set_index("Tissue")["Tissue_Display"].to_dict()
        for tissue, result in tissue_results.items():
            tissue_display = display_lookup[tissue]
            for representation, representation_display in REPRESENTATIONS:
                predictions = result["Predictions"][representation]
                for index, sample_id in enumerate(result["Sample_IDs"]):
                    prediction_rows.append({
                        "Analysis": analysis,
                        "Tissue": tissue,
                        "Tissue_Display": tissue_display,
                        "SAMPID": sample_id,
                        "donor_id": result["Donors"][index],
                        "Fold": result["Folds"][index],
                        "True_SMTSD": result["Labels"][index],
                        "True_Display": result["Class_Display"][result["Labels"][index]],
                        "Representation": representation,
                        "Representation_Display": representation_display,
                        "Predicted_SMTSD": predictions[index],
                        "Predicted_Display": result["Class_Display"][predictions[index]],
                    })
                rows = _complete_confusion_rows(
                    result["Labels"], predictions, result["Classes"], result["Class_Display"],
                )
                for row in rows:
                    confusion_rows.append({
                        "Analysis": analysis,
                        "Tissue": tissue,
                        "Tissue_Display": tissue_display,
                        "Representation": representation,
                        "Representation_Display": representation_display,
                        **row,
                    })

            permutation_rows.extend(_decorate_long_rows(
                result["Permutation_Rows"], analysis, tissue, tissue_display,
            ))
            bootstrap_rows.extend(_decorate_long_rows(
                result["Bootstrap_Rows"], analysis, tissue, tissue_display,
            ))

    predictions_df = pd.DataFrame(prediction_rows)
    assert not predictions_df.duplicated(
        ["Analysis", "Tissue", "SAMPID", "Representation"]
    ).any()
    predictions_df.to_csv(output_dir / "subtissue_oof_predictions.tsv", sep="\t", index=False)
    pd.DataFrame(confusion_rows).to_csv(
        output_dir / "subtissue_confusion_matrices.tsv", sep="\t", index=False,
    )
    pd.DataFrame(permutation_rows).to_csv(
        output_dir / "subtissue_permutation_nulls.tsv", sep="\t", index=False,
    )
    pd.DataFrame(bootstrap_rows).to_csv(
        output_dir / "subtissue_bootstrap_distributions.tsv", sep="\t", index=False,
    )

    count_rows = []
    exclusion_rows = []
    for row in detailed_counts.itertuples(index=False):
        pair = (row.SMTS, row.SMTSD)
        if pair in anatomical_pairs:
            status = "anatomical_primary"
        elif pair in all_pairs:
            status = "supplementary_only"
        elif row.N_Samples < args.min_subtissue_samples:
            status = "excluded_below_min_samples"
        else:
            status = "excluded_no_within_tissue_task"
        count_rows.append({
            "Tissue": row.SMTS,
            "SMTSD": row.SMTSD,
            "Class_Display": class_display_by_tissue[row.SMTS][row.SMTSD],
            "N_Samples": row.N_Samples,
            "N_Donors": row.N_Donors,
            "Included_Anatomical": pair in anatomical_pairs,
            "Included_All_SMTSD": pair in all_pairs,
            "Analysis_Status": status,
        })
        if pair not in anatomical_pairs:
            if row.SMTS == "Blood":
                reason = "transformed_cell_line_vs_primary_tissue_not_anatomical"
            elif row.SMTSD == "Cells - Cultured fibroblasts":
                reason = "cultured_cell_line_not_anatomical"
            elif row.N_Samples < args.min_subtissue_samples:
                reason = "below_min_subtissue_samples"
            else:
                reason = "fewer_than_two_eligible_anatomical_labels"
            exclusion_rows.append({
                "Tissue": row.SMTS,
                "SMTSD": row.SMTSD,
                "N_Samples": row.N_Samples,
                "N_Donors": row.N_Donors,
                "Excluded_From": "anatomical",
                "Reason": reason,
            })
    counts_df = pd.DataFrame(count_rows)
    counts_df.to_csv(output_dir / "subtissue_counts.tsv", sep="\t", index=False)
    pd.DataFrame(exclusion_rows).to_csv(
        output_dir / "subtissue_exclusions.tsv", sep="\t", index=False,
    )

    supplementary_rows = []
    for count_row in count_rows:
        tissue = count_row["Tissue"]
        label = count_row["SMTSD"]
        if (tissue, label) in anatomical_pairs:
            analysis_set = "anatomical"
            result = analysis_task_maps["anatomical"][tissue]
            wide_row = anatomical_results.set_index("Tissue").loc[tissue].to_dict()
        elif (tissue, label) in all_pairs:
            analysis_set = "all_smtsd"
            result = analysis_task_maps["all_smtsd"][tissue]
            wide_row = all_results.set_index("Tissue").loc[tissue].to_dict()
        else:
            analysis_set = "excluded"
            result = None
            wide_row = {}
        supplementary_rows.append({
            **count_row,
            "Analysis_Set": analysis_set,
            "N_Subtissue_Classes": (
                len(result["Classes"]) if result is not None else np.nan
            ),
            "FullLV_Per_Class_Recall": (
                result["Recalls"]["FullLV"][label] if result is not None else np.nan
            ),
            "SHAP_Per_Class_Recall": (
                result["Recalls"]["SHAP"][label] if result is not None else np.nan
            ),
            **{
                key: value
                for key, value in wide_row.items()
                if key.startswith("FullLV_") or key.startswith("SHAP_")
            },
        })
    pd.DataFrame(supplementary_rows).to_csv(
        output_dir / "subtissue_supplementary_table.tsv", sep="\t", index=False,
    )

    dropped_rows = [
        {
            "Tissue": row.SMTS,
            "SMTSD": row.SMTSD,
            "N_Samples": row.N_Samples,
            "Reason": "below_min_subtissue_samples",
        }
        for row in detailed_counts.itertuples(index=False)
        if row.N_Samples < args.min_subtissue_samples
    ]
    _write_legacy_outputs(
        output_dir,
        anatomical_results,
        analysis_task_maps["all_smtsd"],
        rf_accuracy,
        dropped_rows,
    )
    print(f"[gtex] done. results saved to {output_dir}", flush=True)


if __name__ == "__main__":
    main()
