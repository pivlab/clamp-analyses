#!/usr/bin/env python3
# Prepare deterministic four-model drug-repurposing comparison data.
from __future__ import annotations

import argparse
import itertools
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from sklearn.metrics import average_precision_score, roc_auc_score

ORDER = ["ARCHS4", "recount2", "GTEx", "Gene-based"]
METRICS = {"auroc": roc_auc_score, "auprc": average_precision_score}


def bh(pvalues: list[float]) -> np.ndarray:
    values = np.asarray(pvalues, dtype=float)
    order = np.argsort(values)
    ranked = values[order] * len(values) / (np.arange(len(values)) + 1)
    out = np.empty_like(ranked)
    out[order] = np.minimum(np.minimum.accumulate(ranked[::-1])[::-1], 1.0)
    return out


def read_method(path: Path, label: str, gold: pd.DataFrame, n_top: float) -> pd.DataFrame:
    frames = []
    for h5 in sorted(path.glob("*.h5")):
        meta = pd.read_hdf(h5, "metadata")
        if float(meta["n_top_genes"].iloc[0]) != n_top:
            continue
        pred = pd.read_hdf(h5, "prediction").merge(gold, on=["trait", "drug"], how="inner")
        pred["score"] = pred["score"].rank()
        tissue = str(meta["data"].iloc[0]).removeprefix("spredixcan-mashr-zscores-")
        tissue = re.sub(r"-projection-(archs4|recount2|gtex)$|-data$", "", tissue)
        pred["method"], pred["tissue"] = label, tissue
        frames.append(pred[["trait", "drug", "score", "true_class", "method", "tissue"]])
    result = pd.concat(frames, ignore_index=True)
    if result["tissue"].nunique() != 49:
        raise ValueError(f"{label}: expected 49 tissues, got {result['tissue'].nunique()}")
    return result


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--archs4", required=True, type=Path)
    p.add_argument("--recount2", required=True, type=Path)
    p.add_argument("--gtex", required=True, type=Path)
    p.add_argument("--gene", required=True, type=Path)
    p.add_argument("--gold-standard", required=True, type=Path)
    p.add_argument("--n-top-lvs", default=10, type=float)
    p.add_argument("--n-top-genes", default=100, type=float)
    p.add_argument("--n-boot", default=10000, type=int)
    p.add_argument("--seed", default=42, type=int)
    p.add_argument("--output-dir", required=True, type=Path)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    gold = pd.read_pickle(args.gold_standard)
    pred = pd.concat([read_method(args.archs4, "ARCHS4", gold, args.n_top_lvs),
                      read_method(args.recount2, "recount2", gold, args.n_top_lvs),
                      read_method(args.gtex, "GTEx", gold, args.n_top_lvs),
                      read_method(args.gene, "Gene-based", gold, args.n_top_genes)], ignore_index=True)
    per_tissue = pred.groupby(["method", "tissue"], as_index=False).apply(
        lambda x: pd.Series({"n_pairs": len(x), "n_positive": int(x.true_class.sum()),
                             "auroc": roc_auc_score(x.true_class, x.score)}), include_groups=False)
    per_tissue["method"] = pd.Categorical(per_tissue["method"], ORDER, ordered=True)
    per_tissue.sort_values(["method", "tissue"]).to_csv(args.output_dir / "per_tissue_metrics.csv", index=False)
    paired = pred.groupby(["trait", "drug", "method"], as_index=False).agg(score=("score", "max"), true_class=("true_class", "first"))
    wide = paired.pivot(index=["trait", "drug"], columns="method", values="score").reindex(columns=ORDER)
    y = paired.groupby(["trait", "drug"])["true_class"].first().reindex(wide.index).astype(int).to_numpy()
    aggregate = pd.DataFrame({"method": ORDER,
                              **{metric: [func(y, wide[m]) for m in ORDER] for metric, func in METRICS.items()}})
    aggregate.to_csv(args.output_dir / "max_aggregate_reference.csv", index=False)
    rng = np.random.default_rng(args.seed)
    boot = {m: np.empty(args.n_boot) for m in ORDER}
    for i in range(args.n_boot):
        idx = rng.integers(0, len(y), len(y))
        for m in ORDER:
            boot[m][i] = roc_auc_score(y[idx], wide[m].to_numpy()[idx])
    dominance = pd.DataFrame([{"row_m": a, "col_m": b, "value": float(np.mean(boot[a] > boot[b]))}
                              for a, b in itertools.permutations(ORDER, 2)])
    dominance.to_csv(args.output_dir / "ordering_stability.csv", index=False)
    by_method = {m: per_tissue[per_tissue.method == m].set_index("tissue")["auroc"] for m in ORDER}
    tests = []
    for group1, group2 in itertools.combinations(ORDER, 2):
        tests.append({"group1": group1, "group2": group2,
                      "p_value": wilcoxon(by_method[group1], by_method[group2].reindex(by_method[group1].index), alternative="greater").pvalue})
    tests = pd.DataFrame(tests); tests["q_value_bh"] = bh(tests.p_value.tolist())
    tests.to_csv(args.output_dir / "paired_tissue_tests.csv", index=False)
    (args.output_dir / "figure_statistics.txt").write_text(
        f"top_lvs={args.n_top_lvs:g}; top_genes={args.n_top_genes:g}; pairs={len(y)}; positives={int(y.sum())}; tissues=49\n" +
        "\n".join(f"{r.method}: max_AUROC={r.auroc:.4f}, max_AUPRC={r.auprc:.4f}" for r in aggregate.itertuples()) + "\n")


if __name__ == "__main__":
    main()
