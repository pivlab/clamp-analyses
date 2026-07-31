#!/usr/bin/env python3
"""Reconstruct cohort-matched pseudobulk matrices from raw single-cell counts."""

from __future__ import annotations

import argparse
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

from common import (
    h5_gene_names,
    iter_matrix_market,
    load_manifest,
    map_labels,
    matrix_group,
    matrix_market_header,
    obs_values,
    read_csr_rows,
    read_yaml,
    repo_path,
    sparse_shape,
    transform_ids,
)


@dataclass
class Result:
    bulk: pd.DataFrame
    info: pd.DataFrame
    truth_v0: pd.DataFrame
    truth_v1: pd.DataFrame | None
    summary: dict[str, Any]


def count_like(values: np.ndarray) -> dict[str, Any]:
    values = values[values != 0]
    if not len(values):
        return {"nonzero_checked": 0, "count_like_fraction": math.nan}
    return {
        "nonzero_checked": int(len(values)),
        "min_nonzero": float(values.min()),
        "max_nonzero": float(values.max()),
        "count_like_fraction": float(np.isclose(values, np.round(values)).mean()),
    }


def truth_table(samples: np.ndarray, labels: np.ndarray, keep_samples: list[str]) -> pd.DataFrame:
    valid = (samples != "") & (labels != "") & np.isin(samples, keep_samples)
    table = pd.crosstab(samples[valid], labels[valid]).reindex(keep_samples).fillna(0)
    fractions = table.div(table.sum(axis=1).replace(0, np.nan), axis=0).fillna(0)
    fractions.index.name = None
    return fractions


def collapse_duplicate_genes(matrix: np.ndarray, genes: np.ndarray, samples: list[str]) -> pd.DataFrame:
    frame = pd.DataFrame(matrix, index=samples, columns=genes)
    return frame.T.groupby(level=0, sort=True).mean().T


def cohort_selection(
    raw_samples: np.ndarray,
    manifest: pd.DataFrame,
    min_cells: int | None,
) -> tuple[list[str], np.ndarray, pd.Series]:
    requested = manifest["sample"].astype(str).tolist()
    counts = pd.Series(raw_samples).value_counts()
    missing = sorted(set(requested) - set(counts.index))
    if missing:
        raise ValueError(f"{len(missing)} manifest samples are absent from raw metadata: {missing[:10]}")
    if min_cells is not None:
        below = [sample for sample in requested if int(counts[sample]) < min_cells]
        if below:
            raise ValueError(f"manifest samples below min_cells={min_cells}: {below[:10]}")
    mask = np.isin(raw_samples, requested)
    return requested, mask, counts.reindex(requested).astype(int)


def patient_info(manifest: pd.DataFrame, counts: pd.Series) -> pd.DataFrame:
    info = manifest.copy()
    info["sample"] = info["sample"].astype(str)
    info["nCells"] = info["sample"].map(counts).astype(int)
    return info


def aggregate_h5ad(dataset: str, cfg: dict[str, Any], manifest: pd.DataFrame) -> Result:
    started = time.time()
    path = repo_path(cfg["raw"])
    matrix_name = cfg.get("matrix", "X")
    chunk_cells = int(cfg.get("chunk_cells", 5000))
    with h5py.File(path, "r") as h5:
        group = matrix_group(h5, matrix_name)
        n_cells, n_features = sparse_shape(group)
        genes = h5_gene_names(h5, matrix_name)
        if n_features != len(genes):
            raise ValueError(f"matrix has {n_features} features but metadata has {len(genes)}")
        raw_sample = transform_ids(obs_values(h5, cfg["sample_col"]), cfg.get("sample_transform"))
        eligible_sample = raw_sample.copy()
        for column, excluded_values in cfg.get("exclude", {}).items():
            excluded = np.isin(obs_values(h5, column), np.asarray(excluded_values).astype(str))
            eligible_sample[excluded] = ""
        keep_samples, row_mask, sample_counts = cohort_selection(
            eligible_sample, manifest, cfg.get("min_cells")
        )
        sample_codes = pd.Categorical(raw_sample, categories=keep_samples, ordered=True).codes
        raw_v0 = obs_values(h5, cfg["truth_v0_col"])
        labels_v0 = map_labels(raw_v0, cfg.get("truth_v0_map"))
        truth_v0 = truth_table(eligible_sample, labels_v0, keep_samples)
        truth_v1 = None
        if cfg.get("truth_v1_col"):
            raw_v1 = obs_values(h5, cfg["truth_v1_col"])
            labels_v1 = map_labels(raw_v1, cfg.get("truth_v1_map"))
            truth_v1 = truth_table(eligible_sample, labels_v1, keep_samples)

        sums = np.zeros((len(keep_samples), n_features), dtype=np.float64)
        checked: list[np.ndarray] = []
        processed = 0
        for start in range(0, n_cells, chunk_cells):
            stop = min(start + chunk_cells, n_cells)
            keep = row_mask[start:stop]
            if not keep.any():
                continue
            counts = read_csr_rows(group, start, stop, n_features)[keep].tocsr()
            if sum(len(x) for x in checked) < 100_000:
                checked.append(counts.data[: 100_000 - sum(len(x) for x in checked)])
            totals = np.asarray(counts.sum(axis=1)).ravel()
            if np.any(totals <= 0):
                raise ValueError(f"{dataset}: zero-library cells found in retained cohort")
            cpm = counts.multiply(1e6 / totals[:, None]).tocsr()
            codes = sample_codes[start:stop][keep]
            indicator = sp.csr_matrix(
                (np.ones(len(codes)), (codes, np.arange(len(codes)))),
                shape=(len(keep_samples), len(codes)),
            )
            block = (indicator @ cpm).tocoo()
            np.add.at(sums, (block.row, block.col), block.data)
            processed += len(codes)
            if processed % 100_000 < len(codes):
                print(f"[{dataset}] aggregated {processed:,} retained cells", flush=True)

    means = sums / sample_counts.to_numpy()[:, None]
    bulk = collapse_duplicate_genes(means, genes, keep_samples)
    checked_values = np.concatenate(checked) if checked else np.array([])
    summary = {
        "dataset": dataset,
        "raw": str(path),
        "kind": "h5ad",
        "matrix": matrix_name,
        "n_cells_raw": n_cells,
        "n_cells_retained": int(sample_counts.sum()),
        "n_samples_retained": len(keep_samples),
        "n_genes_raw": n_features,
        "n_genes_output": bulk.shape[1],
        "elapsed_seconds": time.time() - started,
        **count_like(checked_values),
    }
    return Result(bulk, patient_info(manifest, sample_counts), truth_v0, truth_v1, summary)


def read_mtx_metadata(cfg: dict[str, Any]) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
    metadata = pd.read_csv(repo_path(cfg["metadata"]), sep="\t", dtype=str)
    barcodes = pd.read_csv(repo_path(cfg["barcodes"]), header=None, dtype=str)[0].to_numpy()
    id_col = cfg["cell_id_col"]
    if not np.array_equal(metadata[id_col].astype(str).to_numpy(), barcodes.astype(str)):
        raise ValueError("Matrix Market metadata and barcode order differ")
    samples = transform_ids(metadata[cfg["sample_col"]], cfg.get("sample_transform"))
    labels_v0 = map_labels(metadata[cfg["truth_v0_col"]].astype(str).to_numpy(), cfg.get("truth_v0_map"))
    labels_v1 = None
    if cfg.get("truth_v1_col"):
        labels_v1 = map_labels(metadata[cfg["truth_v1_col"]].astype(str).to_numpy(), cfg.get("truth_v1_map"))
    return barcodes, samples, labels_v0, labels_v1


def aggregate_mtx(dataset: str, cfg: dict[str, Any], manifest: pd.DataFrame) -> Result:
    started = time.time()
    path = repo_path(cfg["raw"])
    chunk_nnz = int(cfg.get("chunk_nnz", 20_000_000))
    n_genes, n_cells, n_nnz = matrix_market_header(path)
    genes = pd.read_csv(repo_path(cfg["features"]), header=None, dtype=str)[0].to_numpy()
    if n_genes != len(genes):
        raise ValueError("Matrix Market feature count differs from features file")
    _, raw_sample, labels_v0, labels_v1 = read_mtx_metadata(cfg)
    if n_cells != len(raw_sample):
        raise ValueError("Matrix Market cell count differs from metadata")
    keep_samples, cell_mask, sample_counts = cohort_selection(
        raw_sample, manifest, cfg.get("min_cells")
    )
    sample_codes = pd.Categorical(raw_sample, categories=keep_samples, ordered=True).codes
    truth_v0 = truth_table(raw_sample, labels_v0, keep_samples)
    truth_v1 = truth_table(raw_sample, labels_v1, keep_samples) if labels_v1 is not None else None

    totals = np.zeros(n_cells, dtype=np.float64)
    checked: list[np.ndarray] = []
    seen = 0
    for _, col, value in iter_matrix_market(path, chunk_nnz):
        keep = cell_mask[col]
        np.add.at(totals, col[keep], value[keep])
        if sum(len(x) for x in checked) < 100_000:
            checked.append(value[: 100_000 - sum(len(x) for x in checked)])
        seen += len(value)
        print(f"[{dataset}] library-size pass {seen:,}/{n_nnz:,}", flush=True)
    if np.any(totals[cell_mask] <= 0):
        raise ValueError(f"{dataset}: zero-library cells found in retained cohort")

    sums = np.zeros((len(keep_samples), n_genes), dtype=np.float64)
    seen = 0
    for row, col, value in iter_matrix_market(path, chunk_nnz):
        keep = cell_mask[col]
        if keep.any():
            rows = sample_codes[col[keep]]
            values = value[keep] * 1e6 / totals[col[keep]]
            flat = rows.astype(np.int64) * n_genes + row[keep]
            sums += np.bincount(flat, weights=values, minlength=sums.size).reshape(sums.shape)
        seen += len(value)
        print(f"[{dataset}] aggregation pass {seen:,}/{n_nnz:,}", flush=True)

    means = sums / sample_counts.to_numpy()[:, None]
    bulk = collapse_duplicate_genes(means, genes, keep_samples)
    summary = {
        "dataset": dataset,
        "raw": str(path),
        "kind": "matrix_market",
        "n_cells_raw": n_cells,
        "n_cells_retained": int(sample_counts.sum()),
        "n_samples_retained": len(keep_samples),
        "n_genes_raw": n_genes,
        "n_genes_output": bulk.shape[1],
        "n_nonzero_raw": n_nnz,
        "elapsed_seconds": time.time() - started,
        **count_like(np.concatenate(checked)),
    }
    return Result(bulk, patient_info(manifest, sample_counts), truth_v0, truth_v1, summary)


def validate(dataset: str, result: Result) -> None:
    if result.bulk.index.duplicated().any() or result.bulk.columns.duplicated().any():
        raise ValueError(f"{dataset}: duplicate samples or genes")
    values = result.bulk.to_numpy()
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError(f"{dataset}: pseudobulk contains non-finite or negative values")
    if np.any(values.sum(axis=1) <= 0):
        raise ValueError(f"{dataset}: all-zero pseudobulk samples")
    for name, truth in (("truth_v0", result.truth_v0), ("truth_v1", result.truth_v1)):
        if truth is not None and not np.allclose(truth.sum(axis=1), 1):
            raise ValueError(f"{dataset}: {name} rows do not sum to one")
    expected = result.info["sample"].astype(str).tolist()
    if expected != result.bulk.index.astype(str).tolist():
        raise ValueError(f"{dataset}: patient info and expression sample order differ")
    if expected != result.truth_v0.index.astype(str).tolist():
        raise ValueError(f"{dataset}: truth and expression sample order differ")
    if result.summary.get("count_like_fraction", 0) < 0.99:
        raise ValueError(f"{dataset}: selected raw matrix is not count-like")


def write_result(result: Result, args: argparse.Namespace) -> None:
    for path in (args.bulk, args.patient_info, args.truth_v0, args.summary):
        repo_path(path).parent.mkdir(parents=True, exist_ok=True)
    result.bulk.T.to_csv(repo_path(args.bulk))
    result.info.to_csv(repo_path(args.patient_info), index=False)
    result.truth_v0.to_csv(repo_path(args.truth_v0))
    if args.truth_v1:
        if result.truth_v1 is None:
            pd.DataFrame(index=result.bulk.index).to_csv(repo_path(args.truth_v1))
        else:
            result.truth_v1.to_csv(repo_path(args.truth_v1))
    pd.DataFrame([result.summary]).to_csv(repo_path(args.summary), index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--cohort-manifest", required=True)
    parser.add_argument("--bulk", required=True)
    parser.add_argument("--patient-info", required=True)
    parser.add_argument("--truth-v0", required=True)
    parser.add_argument("--truth-v1")
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()
    config = read_yaml(args.config)
    cfg = config["datasets"][args.dataset]
    if not cfg.get("build_pseudobulk", False):
        raise ValueError(
            f"{args.dataset} is configured as a lab-provided pseudobulk input; "
            "only datasets with build_pseudobulk: true may be reconstructed"
        )
    manifest = load_manifest(args.cohort_manifest)
    if cfg["kind"] == "h5ad":
        result = aggregate_h5ad(args.dataset, cfg, manifest)
    elif cfg["kind"] == "matrix_market":
        result = aggregate_mtx(args.dataset, cfg, manifest)
    else:
        raise ValueError(f"unsupported input kind {cfg['kind']!r}")
    validate(args.dataset, result)
    write_result(result, args)
    print(f"[{args.dataset}] completed: {result.bulk.shape[1]} genes x {result.bulk.shape[0]} samples")


if __name__ == "__main__":
    main()
