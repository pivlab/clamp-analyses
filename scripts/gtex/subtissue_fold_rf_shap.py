#!/usr/bin/env python3
# Per-fold donor-grouped RF + out-of-fold SHAP for the subtissue-recovery analysis.
#
# For a single held-out fold: for each broad (SMTS) tissue eligible for subtissue
# recovery, trains a one-vs-rest RandomForest on this fold's training donors
# only (blind to SMTSD subtissue labels throughout), then computes SHAP for the
# held-out fold's tissue-positive samples using the training fold as background.


from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import rpy2.robjects as ro
import shap
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, balanced_accuracy_score
from sklearn.model_selection import RandomizedSearchCV, StratifiedGroupKFold

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import extract_B_matrix  # noqa: E402

readRDS = ro.r["readRDS"]

PARAM_DISTRIBUTIONS = {
    "n_estimators": [100, 300, 500, 1000],
    "max_depth": [None, 10, 20, 30],
    "min_samples_split": [2, 5, 10],
}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--clamp-full-rds", required=True)
    parser.add_argument("--gtex-meta", required=True)
    parser.add_argument("--membership", required=True)
    parser.add_argument("--fold", type=int, required=True)
    parser.add_argument("--tissues", nargs="+", required=True)
    parser.add_argument("--global-seed", type=int, default=42)
    parser.add_argument("--min-rf-accuracy", type=float, default=0.85)
    parser.add_argument("--n-jobs", type=int, required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    np.random.seed(args.global_seed)

    OUTPUT_DIR = Path(args.out_dir)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    gtex_CLAMPfull_rds = readRDS(str(args.clamp_full_rds))
    lv_matrix = extract_B_matrix(gtex_CLAMPfull_rds).T

    # gtex_meta is unused directly here (all label info comes from the fold
    # membership file produced by subtissue_cv_folds.py), but is a declared
    # input for Snakemake staleness tracking, matching sibling scripts.
    _ = pd.read_csv(
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

    membership = pd.read_csv(args.membership, sep="\t")
    membership = membership.set_index("SAMPID")
    lv_matrix = lv_matrix.loc[membership.index]

    fold = args.fold
    train_ids = membership.index[membership["fold"] != fold]
    test_ids = membership.index[membership["fold"] == fold]

    lv_train = lv_matrix.loc[train_ids]
    y_all_train = membership.loc[train_ids, "SMTS"]
    donor_train = membership.loc[train_ids, "donor_id"]

    print(f"[gtex] fold {fold}: {len(train_ids)} train samples, {len(test_ids)} held-out samples", flush=True)

    oof_shap_frames = []
    accuracy_rows = []

    for i, tissue in enumerate(args.tissues):
        print(f"[gtex] fold {fold} [{i + 1}/{len(args.tissues)}] {tissue}", flush=True)

        y_binary_train = (y_all_train == tissue).astype(int)
        n_train_pos = int(y_binary_train.sum())
        n_train_neg = int((1 - y_binary_train).sum())
        imbalance_ratio = n_train_neg / n_train_pos if n_train_pos > 0 else float("inf")

        inner_cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=args.global_seed)
        clf = RandomForestClassifier(random_state=args.global_seed, class_weight="balanced")
        random_search = RandomizedSearchCV(
            estimator=clf,
            param_distributions=PARAM_DISTRIBUTIONS,
            n_iter=15,
            cv=inner_cv,
            scoring="balanced_accuracy",
            n_jobs=args.n_jobs,
            random_state=args.global_seed,
        )
        random_search.fit(lv_train, y_binary_train, groups=donor_train)
        best_params = random_search.best_params_

        final_model = RandomForestClassifier(
            random_state=args.global_seed,
            class_weight="balanced",
            **best_params,
        )
        final_model.fit(lv_train, y_binary_train)

        y_binary_test = (membership.loc[test_ids, "SMTS"] == tissue).astype(int)
        n_test_pos = int(y_binary_test.sum())
        n_test_neg = int((1 - y_binary_test).sum())

        lv_test = lv_matrix.loc[test_ids]
        y_pred_test = final_model.predict(lv_test)
        oof_balanced_accuracy = balanced_accuracy_score(y_binary_test, y_pred_test)
        oof_accuracy = accuracy_score(y_binary_test, y_pred_test)
        passed_gate = bool(oof_balanced_accuracy >= args.min_rf_accuracy)

        accuracy_rows.append({
            "Tissue": tissue,
            "Fold": fold,
            "N_Train_Positive": n_train_pos,
            "N_Train_Negative": n_train_neg,
            "N_Test_Positive": n_test_pos,
            "N_Test_Negative": n_test_neg,
            "Best_Params": str(best_params),
            "OOF_Balanced_Accuracy": oof_balanced_accuracy,
            "OOF_Accuracy": oof_accuracy,
            "Imbalance_Ratio": imbalance_ratio,
            "Passed_Min_Accuracy": passed_gate,
        })
        print(
            f"[gtex]   fold {fold} {tissue}: OOF balanced accuracy = {oof_balanced_accuracy:.4f} "
            f"(passed_gate={passed_gate})",
            flush=True,
        )

        positive_test_ids = test_ids[y_binary_test == 1]
        if len(positive_test_ids) == 0:
            print(f"[gtex]   fold {fold} {tissue}: no held-out positives, skipping SHAP", flush=True)
            continue

        # Background = train fold only (never the held-out fold); explained
        # samples = held-out fold's tissue-positive samples only.
        explainer = shap.TreeExplainer(final_model, data=lv_train, feature_perturbation="interventional")
        shap_explanation = explainer(lv_matrix.loc[positive_test_ids])
        shap_values_class1 = shap_explanation.values[..., 1]

        shap_df = pd.DataFrame(shap_values_class1, index=positive_test_ids, columns=lv_matrix.columns)
        shap_df.index.name = "SAMPID"
        shap_df = shap_df.reset_index()
        shap_df.insert(1, "donor_id", membership.loc[positive_test_ids, "donor_id"].values)
        shap_df.insert(2, "fold", fold)
        shap_df.insert(3, "Tissue", tissue)
        shap_df.insert(4, "SMTSD", membership.loc[positive_test_ids, "SMTSD"].values)
        oof_shap_frames.append(shap_df)

    oof_shap_df = pd.concat(oof_shap_frames, ignore_index=True) if oof_shap_frames else pd.DataFrame()
    accuracy_df = pd.DataFrame(accuracy_rows)

    oof_shap_path = OUTPUT_DIR / "oof_shap.tsv"
    accuracy_path = OUTPUT_DIR / "accuracy_summary.tsv"
    oof_shap_df.to_csv(oof_shap_path, sep="\t", index=False)
    accuracy_df.to_csv(accuracy_path, sep="\t", index=False)
    print(f"[gtex] fold {fold} done. wrote {oof_shap_path}, {accuracy_path}", flush=True)


if __name__ == "__main__":
    main()
