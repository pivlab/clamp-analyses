#!/usr/bin/env python3
"""Project every raw single cell into a full-data CLAMP model without materializing CPM."""

from __future__ import annotations

import argparse
import os
import time
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

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
    sha256_files,
    sparse_shape,
)


def projection_operator(z_path: Path, stats_path: Path, l2_path: Path) -> tuple[pd.Index, np.ndarray, np.ndarray, list[str], float]:
    z_frame = pd.read_csv(z_path, index_col=0)
    stats = pd.read_csv(stats_path, index_col=0)
    genes = z_frame.index.intersection(stats.index, sort=False)
    if len(genes) != len(z_frame):
        missing = z_frame.index.difference(stats.index).tolist()
        raise ValueError(f"row statistics missing {len(missing)} CLAMP genes: {missing[:10]}")
    z = z_frame.loc[genes].to_numpy(dtype=np.float64)
    mean = stats.loc[genes, "mean"].to_numpy(dtype=np.float64)
    variance = stats.loc[genes, "variance"].to_numpy(dtype=np.float64)
    if np.any(~np.isfinite(variance)) or np.any(variance <= 0):
        raise ValueError("CLAMP row statistics contain non-positive variance")
    sigma = np.sqrt(variance)
    l2 = float(pd.read_csv(l2_path)["L2"].iloc[0])
    # Closed-form ridge projection: B = (Y_zscored) @ Z @ (Z'Z + L2*I)^-1, with
    # z-scoring folded into per-gene weights/center so raw counts can be scored
    # directly without ever materializing a z-scored single-cell matrix.
    ridge = np.linalg.inv(z.T @ z + l2 * np.eye(z.shape[1]))
    weights = (z / sigma[:, None]) @ ridge
    center = (mean / sigma) @ z @ ridge
    return genes, weights, center, z_frame.columns.astype(str).tolist(), l2


def raw_feature_weights(raw_genes: np.ndarray, model_genes: pd.Index, model_weights: np.ndarray) -> tuple[np.ndarray, int]:
    # Map each raw feature row to its model weight row; split the weight when
    # a gene ID appears more than once in the raw feature list.
    lookup = {gene: i for i, gene in enumerate(model_genes.astype(str))}
    raw_genes = np.asarray(raw_genes).astype(str)
    counts = pd.Series(raw_genes).value_counts()
    weights = np.zeros((len(raw_genes), model_weights.shape[1]), dtype=np.float64)
    overlap = 0
    for index, gene in enumerate(raw_genes):
        model_index = lookup.get(gene)
        if model_index is not None:
            weights[index] = model_weights[model_index] / int(counts[gene])
            overlap += 1
    if not overlap:
        raise ValueError("raw expression and CLAMP model have no common genes")
    return weights, overlap


def string_dataset(handle: h5py.File, name: str, values: np.ndarray, chunk: int) -> None:
    values = np.asarray(values).astype(str)
    dtype = h5py.string_dtype(encoding="utf-8")
    dataset = handle.create_dataset(name, shape=(len(values),), dtype=dtype, chunks=(min(chunk, len(values)),))
    for start in range(0, len(values), chunk):
        stop = min(start + chunk, len(values))
        dataset[start:stop] = values[start:stop].astype(object)


def create_output(path: Path, n_cells: int, lv_names: list[str], cell_id: np.ndarray,
                  raw_label: np.ndarray, mapped_label: np.ndarray, chunk: int) -> tuple[h5py.File, h5py.Dataset]:
    path.parent.mkdir(parents=True, exist_ok=True)
    handle = h5py.File(path, "w")
    score_chunk = (min(chunk, n_cells), len(lv_names))
    scores = handle.create_dataset(
        "scores", shape=(n_cells, len(lv_names)), dtype="float32", chunks=score_chunk,
        compression="gzip", compression_opts=4, shuffle=True,
    )
    string_dataset(handle, "cell_id", cell_id, chunk)
    string_dataset(handle, "raw_cell_type", raw_label, chunk)
    string_dataset(handle, "mapped_cell_type", mapped_label, chunk)
    string_dataset(handle, "lv_names", np.asarray(lv_names), max(1, len(lv_names)))
    return handle, scores


def label_arrays_h5(h5: h5py.File, cfg: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    cell_id = obs_values(h5)
    raw = obs_values(h5, cfg["raw_cell_type_col"])
    projection_label = obs_values(h5, cfg.get("projection_label_col", cfg["raw_cell_type_col"]))
    mapped = map_labels(projection_label, cfg.get("projection_label_map"))
    return cell_id, raw, mapped


def project_h5ad(dataset: str, cfg: dict, raw_weights: np.ndarray, center: np.ndarray,
                 lv_names: list[str], output: Path, chunk_cells: int) -> dict:
    raw_path = repo_path(cfg["raw"])
    matrix_name = cfg.get("matrix", "X")
    started = time.time()
    checked: list[np.ndarray] = []
    with h5py.File(raw_path, "r") as raw_h5:
        group = matrix_group(raw_h5, matrix_name)
        n_cells, n_features = sparse_shape(group)
        cell_id, raw_label, mapped_label = label_arrays_h5(raw_h5, cfg)
        if len(cell_id) != n_cells or raw_weights.shape[0] != n_features:
            raise ValueError("raw expression dimensions do not match metadata/operator")
        with create_output(output, n_cells, lv_names, cell_id, raw_label, mapped_label, chunk_cells)[0] as out_h5:
            scores = out_h5["scores"]
            # Stream cells in chunks: CPM-normalize each block, then apply the
            # projection operator directly (no full matrix in memory).
            for start in range(0, n_cells, chunk_cells):
                stop = min(start + chunk_cells, n_cells)
                counts = read_csr_rows(group, start, stop, n_features)
                if sum(len(x) for x in checked) < 100_000:
                    checked.append(counts.data[: 100_000 - sum(len(x) for x in checked)])
                totals = np.asarray(counts.sum(axis=1)).ravel()
                if np.any(totals <= 0):
                    raise ValueError(f"{dataset}: zero-library cells in raw matrix")
                block = np.asarray(counts @ raw_weights)
                block *= (1e6 / totals)[:, None]
                block -= center[None, :]
                scores[start:stop] = block.astype(np.float32)
                if start == 0 or stop % (chunk_cells * 10) == 0 or stop == n_cells:
                    print(f"[{dataset}] projected {stop:,}/{n_cells:,} cells", flush=True)
            out_h5.attrs["dataset"] = dataset
            out_h5.attrs["raw_input"] = str(raw_path)
            out_h5.attrs["matrix"] = matrix_name
    values = np.concatenate(checked) if checked else np.array([])
    return {
        "n_cells_raw": n_cells,
        "n_cells_projected": n_cells,
        "n_cells_mapped": int(np.count_nonzero(mapped_label != "")),
        "count_like_fraction": float(np.isclose(values, np.round(values)).mean()),
        "elapsed_seconds": time.time() - started,
        "passes": 1,
    }


def mtx_labels(cfg: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    metadata = pd.read_csv(repo_path(cfg["metadata"]), sep="\t", dtype=str)
    cell_id = metadata[cfg["cell_id_col"]].fillna("").astype(str).to_numpy()
    barcodes = pd.read_csv(repo_path(cfg["barcodes"]), header=None, dtype=str)[0].to_numpy()
    if not np.array_equal(cell_id, barcodes.astype(str)):
        raise ValueError("Matrix Market metadata and barcode order differ")
    raw = metadata[cfg["raw_cell_type_col"]].fillna("").astype(str).to_numpy()
    source = metadata[cfg.get("projection_label_col", cfg["raw_cell_type_col"])].fillna("").astype(str).to_numpy()
    return cell_id, raw, map_labels(source, cfg.get("projection_label_map"))


def project_mtx(dataset: str, cfg: dict, raw_weights: np.ndarray, center: np.ndarray,
                lv_names: list[str], output: Path, chunk_cells: int, chunk_nnz: int) -> dict:
    raw_path = repo_path(cfg["raw"])
    n_features, n_cells, n_nnz = matrix_market_header(raw_path)
    if raw_weights.shape[0] != n_features:
        raise ValueError("Matrix Market feature count does not match projection operator")
    cell_id, raw_label, mapped_label = mtx_labels(cfg)
    if len(cell_id) != n_cells:
        raise ValueError("Matrix Market metadata cell count differs from matrix")
    started = time.time()
    totals = np.zeros(n_cells, dtype=np.float64)
    checked: list[np.ndarray] = []
    seen = 0
    previous_col = -1
    # First pass: accumulate per-cell library sizes (needed for CPM before the
    # per-cell total is known from a single streamed pass).
    for _, col, value in iter_matrix_market(raw_path, chunk_nnz):
        if len(col) and (col[0] < previous_col or np.any(col[1:] < col[:-1])):
            raise ValueError("Matrix Market entries must be cell-column sorted for streaming projection")
        if len(col):
            previous_col = int(col[-1])
        totals += np.bincount(col, weights=value, minlength=n_cells)
        if sum(len(x) for x in checked) < 100_000:
            checked.append(value[: 100_000 - sum(len(x) for x in checked)])
        seen += len(value)
        print(f"[{dataset}] library-size pass {seen:,}/{n_nnz:,}", flush=True)
    if np.any(totals <= 0):
        raise ValueError(f"{dataset}: zero-library cells in raw matrix")

    output.parent.mkdir(parents=True, exist_ok=True)
    memmap_path = output.with_suffix(output.suffix + ".scores.tmp")
    score_map = np.memmap(memmap_path, mode="w+", dtype=np.float64, shape=(n_cells, len(lv_names)))
    for start in range(0, n_cells, chunk_cells):
        score_map[start : start + chunk_cells] = -center
    seen = 0
    previous_col = -1
    # Second pass: apply the projection operator to each streamed NNZ chunk
    for row, col, value in iter_matrix_market(raw_path, chunk_nnz):
        if len(col) and (col[0] < previous_col or np.any(col[1:] < col[:-1])):
            raise ValueError("Matrix Market entries changed order between projection passes")
        if len(col):
            previous_col = int(col[-1])
        scaled = value * 1e6 / totals[col]
        # MTX entries are stored contiguously by cell column. Multiply only the
        # cell span represented by this NNZ chunk; boundary cells are accumulated
        # across chunks. This is the same global operator, with bounded RAM.
        first, stop = int(col[0]), int(col[-1]) + 1
        counts = sp.coo_matrix(
            (scaled, (col - first, row)), shape=(stop - first, n_features)
        ).tocsr()
        score_map[first:stop] += np.asarray(counts @ raw_weights)
        score_map.flush()
        seen += len(value)
        print(f"[{dataset}] projection pass {seen:,}/{n_nnz:,}", flush=True)

    handle, scores = create_output(output, n_cells, lv_names, cell_id, raw_label, mapped_label, chunk_cells)
    with handle:
        for start in range(0, n_cells, chunk_cells):
            stop = min(start + chunk_cells, n_cells)
            scores[start:stop] = score_map[start:stop].astype(np.float32)
        handle.attrs["dataset"] = dataset
        handle.attrs["raw_input"] = str(raw_path)
        handle.attrs["matrix"] = "matrix_market"
    del score_map
    os.unlink(memmap_path)
    values = np.concatenate(checked)
    return {
        "n_cells_raw": n_cells,
        "n_cells_projected": n_cells,
        "n_cells_mapped": int(np.count_nonzero(mapped_label != "")),
        "count_like_fraction": float(np.isclose(values, np.round(values)).mean()),
        "elapsed_seconds": time.time() - started,
        "passes": 2,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--z", required=True)
    parser.add_argument("--l2", required=True)
    parser.add_argument("--row-stats", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--chunk-cells", type=int, default=20_000)
    parser.add_argument("--chunk-nnz", type=int, default=20_000_000)
    args = parser.parse_args()
    config = read_yaml(args.config)
    cfg = config["datasets"][args.dataset]
    model_genes, model_weights, center, lv_names, l2 = projection_operator(
        repo_path(args.z), repo_path(args.row_stats), repo_path(args.l2)
    )
    if cfg["kind"] == "h5ad":
        with h5py.File(repo_path(cfg["raw"]), "r") as raw_h5:
            raw_genes = h5_gene_names(raw_h5, cfg.get("matrix", "X"))
    else:
        raw_genes = pd.read_csv(repo_path(cfg["features"]), header=None, dtype=str)[0].to_numpy()
    feature_weights, overlap = raw_feature_weights(raw_genes, model_genes, model_weights)
    output = repo_path(args.output)
    if cfg["kind"] == "h5ad":
        summary = project_h5ad(args.dataset, cfg, feature_weights, center, lv_names, output, args.chunk_cells)
    elif cfg["kind"] == "matrix_market":
        summary = project_mtx(
            args.dataset, cfg, feature_weights, center, lv_names, output,
            args.chunk_cells, args.chunk_nnz,
        )
    else:
        raise ValueError(f"unsupported input kind {cfg['kind']!r}")
    summary.update({
        "dataset": args.dataset,
        "n_model_genes": len(model_genes),
        "n_raw_features_overlapping": overlap,
        "n_lvs": len(lv_names),
        "L2": l2,
        "normalization": "full_pseudobulk_training_statistics",
        "model_preprocess_sha256": sha256_files(args.z, args.l2, args.row_stats),
    })
    with h5py.File(output, "r+") as out_h5:
        out_h5.attrs["model"] = "CLAMPfull"
        out_h5.attrs["model_z"] = str(repo_path(args.z))
        out_h5.attrs["model_l2"] = str(repo_path(args.l2))
        out_h5.attrs["preprocessing_row_stats"] = str(repo_path(args.row_stats))
        out_h5.attrs["model_preprocess_sha256"] = summary["model_preprocess_sha256"]
        out_h5.attrs["normalization"] = summary["normalization"]
        out_h5.attrs["n_model_genes"] = summary["n_model_genes"]
        out_h5.attrs["n_raw_features_overlapping"] = summary["n_raw_features_overlapping"]
        out_h5.attrs["n_cells_projected"] = summary["n_cells_projected"]
        out_h5.attrs["L2"] = summary["L2"]
        out_h5.attrs["chunk_cells"] = args.chunk_cells
        out_h5.attrs["chunk_nnz"] = args.chunk_nnz
        out_h5.attrs["completed_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    if summary["count_like_fraction"] < 0.99:
        raise ValueError(f"{args.dataset}: raw expression is not count-like")
    summary_path = repo_path(args.summary)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([summary]).to_csv(summary_path, index=False)
    print(f"[{args.dataset}] wrote {summary['n_cells_projected']:,} cell projections to {output}")


if __name__ == "__main__":
    main()
