#!/usr/bin/env python3
"""Build donor-bulk libraries and matched sampled single-cell matrices."""

from __future__ import annotations

import argparse
import sys
import tempfile
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts" / "pseudobulk"))

from common import (
    h5_gene_names,
    iter_matrix_market,
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


def validate_values(dataset: str, values: np.ndarray) -> None:
    if (
        not np.isfinite(values).all()
        or np.any(values < 0)
        or not np.allclose(values, np.round(values))
    ):
        raise ValueError(f"{dataset}: raw expression is not finite, nonnegative integer-like data")


def load_manifest(cfg: dict, manifest_path: str | None = None) -> pd.DataFrame:
    if manifest_path:
        candidates = [repo_path(manifest_path)]
    else:
        candidates = []
    if cfg.get("cohort_manifest"):
        candidates.append(repo_path(cfg["cohort_manifest"]))
    dataset = cfg["dataset"]
    candidates.extend(
        [
            repo_path(f"output/01_model_building/00_pseudobulk/{dataset}/pseudobulk/cohort_manifest.csv"),
            repo_path(f"output/01_model_building/00_pseudobulk/{dataset}/pseudobulk/patient_info.csv"),
            repo_path(f"data/pseudobulk/{dataset}/patient_info.csv"),
        ]
    )
    for path in candidates:
        if path.exists():
            frame = pd.read_csv(path, dtype={"sample": str})
            if "sample" not in frame or frame["sample"].duplicated().any():
                raise ValueError(f"malformed cohort manifest: {path}")
            return frame
    raise FileNotFoundError(f"{dataset}: no cohort manifest found")


def select_cells(
    labels: np.ndarray,
    eligible: np.ndarray,
    target: int,
    minimum_per_type: int,
    seed: int,
) -> np.ndarray:
    retained = np.flatnonzero(eligible)
    if not len(retained):
        raise ValueError("no retained cells available for expression UMAP")
    if len(retained) <= target:
        return retained.astype(np.int64, copy=False)
    rng = np.random.default_rng(seed)
    selected: set[int] = set()
    for cell_type in pd.unique(labels[eligible]):
        observed = np.flatnonzero(eligible & (labels == cell_type))
        quota = min(minimum_per_type, len(observed))
        selected.update(rng.choice(observed, size=quota, replace=False).tolist())
    if len(selected) > target:
        raise ValueError(
            f"minimum-per-type quotas require {len(selected)} cells, above target {target}"
        )
    remaining = target - len(selected)
    if remaining:
        pool = np.setdiff1d(retained, np.fromiter(selected, dtype=np.int64), assume_unique=False)
        selected.update(rng.choice(pool, size=remaining, replace=False).tolist())
    result = np.asarray(sorted(selected), dtype=np.int64)
    if len(result) != target:
        raise ValueError("fixed-size stratified sampling failed")
    return result


def donor_codes(
    raw_samples: np.ndarray,
    manifest: pd.DataFrame,
    eligible: np.ndarray,
) -> tuple[list[str], np.ndarray, np.ndarray, np.ndarray]:
    samples = manifest["sample"].astype(str).tolist()
    missing = sorted(set(samples) - set(raw_samples[eligible]))
    if missing:
        raise ValueError(f"manifest donors absent after cell filtering: {missing[:10]}")
    codes = pd.Categorical(raw_samples, categories=samples, ordered=True).codes
    keep = eligible & (codes >= 0)
    counts = np.bincount(codes[keep], minlength=len(samples)).astype(np.int64)
    if np.any(counts <= 0):
        raise ValueError("at least one retained donor contains no eligible cells")
    return samples, codes, keep, counts


def truth_table(
    samples: np.ndarray,
    labels: np.ndarray,
    keep: np.ndarray,
    donor_order: list[str],
) -> pd.DataFrame:
    table = pd.crosstab(samples[keep], labels[keep]).reindex(donor_order).fillna(0)
    truth = table.div(table.sum(axis=1), axis=0)
    truth.index.name = None
    return truth


def gene_aggregator(genes: np.ndarray) -> tuple[np.ndarray, sp.csr_matrix]:
    genes = np.asarray(genes).astype(str)
    unique, codes = np.unique(genes, return_inverse=True)
    aggregator = sp.csr_matrix(
        (np.ones(len(genes), dtype=np.float64), (np.arange(len(genes)), codes)),
        shape=(len(genes), len(unique)),
    )
    return unique, aggregator


def metadata_arrays_h5(h5: h5py.File, cfg: dict) -> tuple[np.ndarray, ...]:
    n_cells = len(obs_values(h5))
    eligible = np.ones(n_cells, dtype=bool)
    for column, excluded in cfg.get("exclude", {}).items():
        eligible &= ~np.isin(obs_values(h5, column), np.asarray(excluded).astype(str))
    raw_samples = transform_ids(obs_values(h5, cfg["sample_col"]), cfg.get("sample_transform"))
    raw_labels = obs_values(h5, cfg["truth_v0_col"])
    labels = map_labels(raw_labels, cfg.get("truth_v0_map"))
    eligible &= labels != ""
    return obs_values(h5), raw_samples, labels, eligible


def aggregate_h5(
    dataset: str,
    cfg: dict,
    manifest: pd.DataFrame,
    chunk_cells: int,
    target_cells: int,
    minimum_per_type: int,
    seed: int,
) -> tuple:
    raw_path = repo_path(cfg["raw"])
    matrix_name = cfg.get("matrix", "X")
    with h5py.File(raw_path, "r") as h5:
        group = matrix_group(h5, matrix_name)
        n_cells, n_genes = sparse_shape(group)
        n_nnz = int(group["data"].shape[0])
        genes = h5_gene_names(h5, matrix_name)
        cell_ids, raw_samples, labels, eligible = metadata_arrays_h5(h5, cfg)
        if len(cell_ids) != n_cells:
            raise ValueError(f"{dataset}: expression and metadata cell counts differ")
        donors, codes, keep, cell_counts = donor_codes(raw_samples, manifest, eligible)
        truth = truth_table(raw_samples, labels, keep, donors)
        sample_index = select_cells(labels, keep, target_cells, minimum_per_type, seed)
        selected_lookup = np.zeros(n_cells, dtype=bool)
        selected_lookup[sample_index] = True
        sums = np.zeros((len(donors), n_genes), dtype=np.float64)
        sampled_blocks = []
        retained_processed = 0
        for start in range(0, n_cells, chunk_cells):
            stop = min(start + chunk_cells, n_cells)
            block = read_csr_rows(group, start, stop, n_genes).tocsr()
            validate_values(dataset, block.data)
            local_sample = selected_lookup[start:stop]
            if local_sample.any():
                sampled_blocks.append(block[local_sample])
            local_keep = keep[start:stop]
            if local_keep.any():
                retained = block[local_keep]
                totals = np.asarray(retained.sum(axis=1)).ravel()
                if np.any(totals <= 0):
                    raise ValueError(f"{dataset}: retained zero-library cell")
                local_codes = codes[start:stop][local_keep]
                indicator = sp.csr_matrix(
                    (np.ones(len(local_codes)), (local_codes, np.arange(len(local_codes)))),
                    shape=(len(donors), len(local_codes)),
                )
                donor_block = (indicator @ retained).tocoo()
                np.add.at(sums, (donor_block.row, donor_block.col), donor_block.data)
                retained_processed += int(local_keep.sum())
            if stop % (chunk_cells * 10) == 0 or stop == n_cells:
                print(f"[{dataset}] checked {stop:,}/{n_cells:,} cells", flush=True)
        sampled = sp.vstack(sampled_blocks, format="csr")
    summary = {
        "dataset": dataset,
        "kind": "h5ad",
        "raw": str(raw_path),
        "matrix": matrix_name,
        "n_cells_raw": n_cells,
        "n_cells_retained": int(cell_counts.sum()),
        "n_cells_sampled": len(sample_index),
        "n_donors": len(donors),
        "n_genes_raw": n_genes,
        "n_nonzero_raw": n_nnz,
        "raw_total_counts_retained": float(sums.sum()),
        "raw_validation_scope": "exhaustive_all_nonzeros",
    }
    return (
        sums,
        genes,
        truth,
        donors,
        cell_counts,
        sampled,
        sample_index,
        cell_ids[sample_index],
        raw_samples[sample_index],
        labels[sample_index],
        summary,
    )


def metadata_arrays_mtx(cfg: dict, n_cells: int) -> tuple[np.ndarray, ...]:
    metadata = pd.read_csv(repo_path(cfg["metadata"]), sep="\t", dtype=str).fillna("")
    cell_ids = metadata[cfg["cell_id_col"]].astype(str).to_numpy()
    barcodes = pd.read_csv(repo_path(cfg["barcodes"]), header=None, dtype=str)[0].to_numpy()
    if len(cell_ids) != n_cells or not np.array_equal(cell_ids, barcodes.astype(str)):
        raise ValueError("Matrix Market metadata and barcode order differ")
    eligible = np.ones(n_cells, dtype=bool)
    for column, excluded in cfg.get("exclude", {}).items():
        eligible &= ~metadata[column].astype(str).isin(np.asarray(excluded).astype(str)).to_numpy()
    raw_samples = transform_ids(metadata[cfg["sample_col"]], cfg.get("sample_transform"))
    labels = map_labels(metadata[cfg["truth_v0_col"]].astype(str).to_numpy(), cfg.get("truth_v0_map"))
    eligible &= labels != ""
    return cell_ids, raw_samples, labels, eligible


def aggregate_mtx(
    dataset: str,
    cfg: dict,
    manifest: pd.DataFrame,
    chunk_nnz: int,
    target_cells: int,
    minimum_per_type: int,
    seed: int,
    temp_parent: Path,
) -> tuple:
    raw_path = repo_path(cfg["raw"])
    n_genes, n_cells, n_nnz = matrix_market_header(raw_path)
    genes = pd.read_csv(repo_path(cfg["features"]), header=None, dtype=str)[0].to_numpy()
    cell_ids, raw_samples, labels, eligible = metadata_arrays_mtx(cfg, n_cells)
    donors, codes, keep, cell_counts = donor_codes(raw_samples, manifest, eligible)
    truth = truth_table(raw_samples, labels, keep, donors)
    sample_index = select_cells(labels, keep, target_cells, minimum_per_type, seed)
    sample_lookup = np.full(n_cells, -1, dtype=np.int32)
    sample_lookup[sample_index] = np.arange(len(sample_index), dtype=np.int32)
    totals = np.zeros(n_cells, dtype=np.float64)
    sample_nnz = 0
    seen = 0
    for _, col, value in iter_matrix_market(raw_path, chunk_nnz):
        validate_values(dataset, value)
        np.add.at(totals, col, value)
        sample_nnz += int(np.count_nonzero(sample_lookup[col] >= 0))
        seen += len(value)
        print(f"[{dataset}] validation pass {seen:,}/{n_nnz:,}", flush=True)
    if np.any(totals[keep] <= 0):
        raise ValueError(f"{dataset}: retained zero-library cell")

    sums = np.zeros((len(donors), n_genes), dtype=np.float64)
    temp_parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=f"{dataset}_sample_", dir=temp_parent) as temp_dir:
        temp = Path(temp_dir)
        sample_rows = np.memmap(temp / "rows.bin", mode="w+", dtype=np.int32, shape=(sample_nnz,))
        sample_cols = np.memmap(temp / "cols.bin", mode="w+", dtype=np.int32, shape=(sample_nnz,))
        sample_data = np.memmap(temp / "data.bin", mode="w+", dtype=np.float64, shape=(sample_nnz,))
        cursor = 0
        seen = 0
        for row, col, value in iter_matrix_market(raw_path, chunk_nnz):
            retained = keep[col]
            if retained.any():
                kept_col = col[retained]
                flat = codes[kept_col].astype(np.int64) * n_genes + row[retained]
                sums += np.bincount(flat, weights=value[retained], minlength=sums.size).reshape(sums.shape)
            local_rows = sample_lookup[col]
            selected = local_rows >= 0
            count = int(selected.sum())
            if count:
                sample_rows[cursor : cursor + count] = local_rows[selected]
                sample_cols[cursor : cursor + count] = row[selected]
                sample_data[cursor : cursor + count] = value[selected]
                cursor += count
            seen += len(value)
            print(f"[{dataset}] aggregation pass {seen:,}/{n_nnz:,}", flush=True)
        if cursor != sample_nnz:
            raise ValueError(f"{dataset}: sampled NNZ accounting mismatch")
        sampled = sp.coo_matrix(
            (sample_data, (sample_rows, sample_cols)),
            shape=(len(sample_index), n_genes),
        ).tocsr()
        del sample_rows, sample_cols, sample_data
    summary = {
        "dataset": dataset,
        "kind": "matrix_market",
        "raw": str(raw_path),
        "matrix": "matrix_market",
        "n_cells_raw": n_cells,
        "n_cells_retained": int(cell_counts.sum()),
        "n_cells_sampled": len(sample_index),
        "n_donors": len(donors),
        "n_genes_raw": n_genes,
        "n_nonzero_raw": n_nnz,
        "raw_total_counts_retained": float(sums.sum()),
        "raw_validation_scope": "exhaustive_all_nonzeros",
    }
    return (
        sums,
        genes,
        truth,
        donors,
        cell_counts,
        sampled,
        sample_index,
        cell_ids[sample_index],
        raw_samples[sample_index],
        labels[sample_index],
        summary,
    )


def write_outputs(args: argparse.Namespace, result: tuple, manifest: pd.DataFrame) -> None:
    (
        sums,
        genes,
        truth,
        donors,
        cell_counts,
        sampled,
        sample_index,
        sample_ids,
        sample_donors,
        sample_labels,
        summary,
    ) = result
    unique_genes, aggregator = gene_aggregator(genes)
    collapsed = np.asarray(sums @ aggregator)
    collapsed_sample = (sampled @ aggregator).tocsr()
    if not np.isclose(collapsed.sum(), sums.sum()):
        raise ValueError("duplicate-gene collapse changed total donor counts")
    if not np.isclose(collapsed_sample.sum(), sampled.sum()):
        raise ValueError("duplicate-gene collapse changed sampled-cell counts")
    if list(truth.index.astype(str)) != donors or not np.allclose(truth.sum(axis=1), 1):
        raise ValueError("truth rows and donor matrix columns are not aligned")
    implied = truth.to_numpy() * cell_counts[:, None]
    if (
        not np.allclose(implied, np.round(implied), atol=1e-7)
        or not np.array_equal(np.round(implied).sum(axis=1).astype(np.int64), cell_counts)
    ):
        raise ValueError("truth fractions do not reproduce retained-cell denominators")

    paths = [
        args.counts,
        args.truth,
        args.patient_info,
        args.summary,
        args.sample_counts,
        args.sample_cells,
        args.sample_genes,
    ]
    for path in paths:
        repo_path(path).parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(collapsed.T, index=unique_genes, columns=donors).to_csv(repo_path(args.counts))
    truth.to_csv(repo_path(args.truth))
    info = manifest.set_index(manifest["sample"].astype(str)).loc[donors].reset_index(drop=True).copy()
    info["sample"] = donors
    info["nCells"] = cell_counts
    info.to_csv(repo_path(args.patient_info), index=False)
    sp.save_npz(repo_path(args.sample_counts), collapsed_sample, compressed=True)
    pd.DataFrame({"gene": unique_genes}).to_csv(repo_path(args.sample_genes), index=False)
    pd.DataFrame(
        {
            "sample_row": np.arange(len(sample_index)),
            "cell_index": sample_index,
            "cell_id": sample_ids,
            "sample": sample_donors,
            "mapped_cell_type": sample_labels,
        }
    ).to_csv(repo_path(args.sample_cells), index=False)
    summary.update(
        {
            "n_genes_output": len(unique_genes),
            "n_nonzero_sampled": int(collapsed_sample.nnz),
            "collapsed_total_counts": float(collapsed.sum()),
            "truth_denominators_equal_expression_cells": True,
            "min_donor_library": float(collapsed.sum(axis=1).min()),
            "max_donor_library": float(collapsed.sum(axis=1).max()),
            "sample_seed": args.seed,
            "sample_target": args.sample_cells_target,
            "sample_minimum_per_type": args.minimum_per_type,
        }
    )
    pd.DataFrame([summary]).to_csv(repo_path(args.summary), index=False)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--counts", required=True)
    parser.add_argument("--truth", required=True)
    parser.add_argument("--patient-info", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--sample-counts", required=True)
    parser.add_argument("--sample-cells", required=True)
    parser.add_argument("--sample-genes", required=True)
    parser.add_argument("--chunk-cells", type=int, default=20_000)
    parser.add_argument("--chunk-nnz", type=int, default=20_000_000)
    parser.add_argument("--sample-cells-target", type=int, default=60_000)
    parser.add_argument("--minimum-per-type", type=int, default=250)
    parser.add_argument("--seed", type=int, default=123)
    args = parser.parse_args()
    config = read_yaml(args.config)
    cfg = dict(config["datasets"][args.dataset], dataset=args.dataset)
    manifest = load_manifest(cfg, args.manifest)
    started = time.time()
    if cfg["kind"] == "h5ad":
        result = aggregate_h5(
            args.dataset, cfg, manifest, args.chunk_cells,
            args.sample_cells_target, args.minimum_per_type, args.seed,
        )
    elif cfg["kind"] == "matrix_market":
        result = aggregate_mtx(
            args.dataset, cfg, manifest, args.chunk_nnz,
            args.sample_cells_target, args.minimum_per_type, args.seed,
            repo_path(args.sample_counts).parent,
        )
    else:
        raise ValueError(f"unsupported input kind: {cfg['kind']}")
    result[-1]["elapsed_seconds"] = time.time() - started
    write_outputs(args, result, manifest)
    print(f"[{args.dataset}] donor-bulk aggregation complete", flush=True)


if __name__ == "__main__":
    main()
