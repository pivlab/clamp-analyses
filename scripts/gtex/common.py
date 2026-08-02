"""Shared helpers for GTEx production scripts."""

from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import balanced_accuracy_score
from sklearn.model_selection import RandomizedSearchCV, StratifiedKFold


def derive_donor_id(sampid: pd.Series) -> pd.Series:
    """Derive the GTEx donor/subject ID from SAMPID (first two hyphen-delimited
    components, e.g. 'GTEX-1117F-0226-SM-5GZZ7' -> 'GTEX-1117F'). Matches SUBJID
    in GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt.
    """
    return sampid.astype(str).str.split("-").str[:2].str.join("-")


def extract_B_matrix(rds_obj) -> pd.DataFrame:
    """Extract a B matrix from an RDS object (readRDS result) into a DataFrame."""
    B_matrix = rds_obj.rx2("B")
    with localconverter(ro.default_converter + pandas2ri.converter):
        B_values = ro.conversion.rpy2py(B_matrix)
    return pd.DataFrame(
        data=B_values,
        index=B_matrix.rownames if B_matrix.rownames else None,
        columns=B_matrix.colnames if B_matrix.colnames else None,
    )


def extract_summary_matrix(rds_obj) -> pd.DataFrame:
    """Extract a CLAMPfull summary table from an RDS object into a DataFrame."""
    summary_matrix = rds_obj.rx2("summary")
    with localconverter(ro.default_converter + pandas2ri.converter):
        summary_values = ro.conversion.rpy2py(summary_matrix)
    return pd.DataFrame(
        data=summary_values,
        index=summary_matrix.rownames if summary_matrix.rownames else None,
        columns=summary_matrix.colnames if summary_matrix.colnames else None,
    )


def run_tissue_rf_shap(
    lv_matrix: pd.DataFrame,
    tissues: list[str],
    label_fn: Callable[[str], pd.Series],
    *,
    output_dir: Path,
    global_seed: int,
    min_rf_accuracy: float,
    param_distributions: dict,
    extra_accuracy_fields_fn: Callable[[str], dict] | None = None,
) -> None:
    """Per-tissue nested-CV RandomForest + SHAP training loop.

    Body copied verbatim from 02_rf_kmeans/00_LV_importance_kmeans.ipynb and
    03_rf_true_labels/00_LV_importance_true_labels.ipynb, which were
    otherwise near-identical: label_fn and extra_accuracy_fields_fn are the
    only seams, letting each caller supply its own tissue-label derivation
    (k-means cluster majority vote vs. raw GTEx SMTSD) and its own extra
    accuracy-row columns.
    """
    (output_dir / "per_tissue").mkdir(parents=True, exist_ok=True)

    accuracy_list = []
    shap_results_list = []
    cumulative_importance_list = []
    top5_lvs_by_tissue = {}

    for i, tissue in enumerate(tissues):
        print(f"[gtex] [{i + 1}/{len(tissues)}] {tissue}", flush=True)

        y_binary = label_fn(tissue).astype(int)
        n_pos = y_binary.sum()
        n_neg = (1 - y_binary).sum()
        imbalance_ratio = n_neg / n_pos if n_pos > 0 else float("inf")

        outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=global_seed + i)
        outer_test_scores = []
        best_params_list = []

        for outer_fold, (dev_idx, test_idx) in enumerate(outer_cv.split(lv_matrix, y_binary)):
            lv_data_dev = lv_matrix.iloc[dev_idx]
            y_dev = y_binary.iloc[dev_idx]
            lv_data_test = lv_matrix.iloc[test_idx]
            y_test = y_binary.iloc[test_idx]

            inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=global_seed + i)
            clf = RandomForestClassifier(random_state=global_seed + i, class_weight="balanced")
            random_search = RandomizedSearchCV(
                estimator=clf,
                param_distributions=param_distributions,
                n_iter=15,
                cv=inner_cv,
                scoring="balanced_accuracy",
                n_jobs=-1,
                random_state=global_seed + i,
            )
            random_search.fit(lv_data_dev, y_dev)
            best_params_list.append(random_search.best_params_)

            final_model_outer = RandomForestClassifier(
                random_state=global_seed + i,
                class_weight="balanced",
                **random_search.best_params_,
            )
            final_model_outer.fit(lv_data_dev, y_dev)
            y_pred = final_model_outer.predict(lv_data_test)
            outer_test_scores.append(balanced_accuracy_score(y_test, y_pred))

        best_fold_idx = int(np.argmax(outer_test_scores))
        final_params = best_params_list[best_fold_idx]

        print(
            f"[gtex]   nested CV balanced accuracy: "
            f"{np.mean(outer_test_scores):.4f} +/- {np.std(outer_test_scores):.4f}",
            flush=True,
        )

        final_model = RandomForestClassifier(
            random_state=global_seed + i,
            class_weight="balanced",
            **final_params,
        )
        final_model.fit(lv_matrix, y_binary)

        import shap

        explainer = shap.TreeExplainer(final_model, data=lv_matrix, feature_perturbation="interventional")
        shap_explanation = explainer(lv_matrix)
        shap_values_class1 = shap_explanation.values[..., 1]

        tissue_mask = y_binary == 1
        shap_values_tissue_samples = shap_values_class1[tissue_mask]
        mean_shap_tissue = np.mean(shap_values_tissue_samples, axis=0)

        df_shap_all = pd.DataFrame({
            "Feature": lv_matrix.columns,
            "Mean_SHAP_Tissue": mean_shap_tissue,
        })
        df_positive = df_shap_all[df_shap_all["Mean_SHAP_Tissue"] > 0].copy()
        df_positive = df_positive.sort_values("Mean_SHAP_Tissue", ascending=False).reset_index(drop=True)

        n_features_needed = {}
        if len(df_positive) > 0:
            total_positive_shap = df_positive["Mean_SHAP_Tissue"].sum()
            df_positive["Cumulative_SHAP"] = df_positive["Mean_SHAP_Tissue"].cumsum()
            df_positive["Cumulative_Percent"] = (df_positive["Cumulative_SHAP"] / total_positive_shap) * 100
            df_positive["Rank"] = range(1, len(df_positive) + 1)

            for thresh in [50, 70, 80, 90, 95]:
                n_features_needed[thresh] = (df_positive["Cumulative_Percent"] >= thresh).idxmax() + 1

            df_cumulative = df_positive[["Feature", "Mean_SHAP_Tissue", "Cumulative_Percent", "Rank"]].copy()
            df_cumulative["Tissue"] = tissue
            cumulative_importance_list.append(df_cumulative)

        accuracy_row = {
            "Tissue": tissue,
            **(extra_accuracy_fields_fn(tissue) if extra_accuracy_fields_fn else {}),
            "Best_Test_Balanced_Accuracy": outer_test_scores[best_fold_idx],
            "Mean_CV_Accuracy": np.mean(outer_test_scores),
            "Std_CV_Accuracy": np.std(outer_test_scores),
            "Imbalance_Ratio": imbalance_ratio,
            "N_Positive_LVs": len(df_positive),
            "LVs_for_80pct": n_features_needed.get(80),
            "LVs_for_90pct": n_features_needed.get(90),
        }
        accuracy_list.append(accuracy_row)

        df_all_positive = df_positive.copy()
        df_all_positive["Tissue"] = tissue
        shap_results_list.append(df_all_positive)

        safe_name = tissue.replace(" ", "_").replace("/", "_")
        df_positive.to_csv(output_dir / "per_tissue" / f"shap_positive_{safe_name}.tsv", sep="\t", index=False)
        top5_lvs_by_tissue[tissue] = df_positive["Feature"].head(5).tolist()
        print(f"[gtex]   top LVs: {top5_lvs_by_tissue[tissue]}", flush=True)

    accuracy_df = pd.DataFrame(accuracy_list)
    shap_results_df = pd.concat(shap_results_list, ignore_index=True)
    cumulative_df = pd.concat(cumulative_importance_list, ignore_index=True)

    passing_tissues = accuracy_df[accuracy_df["Best_Test_Balanced_Accuracy"] >= min_rf_accuracy]["Tissue"]
    accuracy_df = accuracy_df[accuracy_df["Tissue"].isin(passing_tissues)]
    shap_results_df = shap_results_df[shap_results_df["Tissue"].isin(passing_tissues)]
    cumulative_df = cumulative_df[cumulative_df["Tissue"].isin(passing_tissues)]
    print(f"[gtex] Tissues passing accuracy >= {min_rf_accuracy}: {sorted(passing_tissues.tolist())}", flush=True)

    accuracy_df.to_csv(output_dir / "accuracy_summary.tsv", sep="\t", index=False)
    shap_results_df.to_csv(output_dir / "all_shap_positive.tsv", sep="\t", index=False)
    cumulative_df.to_csv(output_dir / "cumulative_importance.tsv", sep="\t", index=False)
    print(f"[gtex] done. results saved to: {output_dir}", flush=True)
