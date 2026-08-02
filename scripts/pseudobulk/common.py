"""Shared low-level helpers for pseudobulk production scripts."""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path
from typing import Any, Iterator

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import yaml


ROOT = Path(__file__).resolve().parents[2]


def repo_path(path: str | Path) -> Path:
    value = Path(path)
    return value if value.is_absolute() else ROOT / value


def read_yaml(path: str | Path) -> dict[str, Any]:
    with open(repo_path(path), encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def sha256_files(*paths: str | Path) -> str:
    digest = hashlib.sha256()
    for path in paths:
        with open(repo_path(path), "rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
    return digest.hexdigest()


def decode_array(values: np.ndarray) -> np.ndarray:
    if values.dtype.kind in {"S", "O"}:
        return np.asarray(
            [x.decode() if isinstance(x, bytes) else str(x) for x in values],
            dtype=str,
        )
    return values.astype(str)


def read_h5_string(node: h5py.Dataset | h5py.Group) -> np.ndarray:
    if isinstance(node, h5py.Group):
        codes = node["codes"][:]
        categories = decode_array(node["categories"][:])
        result = np.full(len(codes), "", dtype=object)
        valid = (codes >= 0) & (codes < len(categories))
        result[valid] = categories[codes[valid]]
        return result.astype(str)
    return decode_array(node[:])


def obs_values(h5: h5py.File, column: str | None = None) -> np.ndarray:
    obs = h5["obs"]
    key = column or obs.attrs.get("_index", "_index")
    if isinstance(key, bytes):
        key = key.decode()
    if key not in obs:
        raise KeyError(f"obs column {key!r} not found in {h5.filename}")
    return read_h5_string(obs[key])


def matrix_group(h5: h5py.File, matrix_name: str) -> h5py.Group:
    key = "raw/X" if matrix_name == "raw.X" else matrix_name
    if key not in h5:
        raise KeyError(f"matrix {matrix_name!r} not found in {h5.filename}")
    group = h5[key]
    if not isinstance(group, h5py.Group):
        raise TypeError(f"matrix {matrix_name!r} is not a backed sparse matrix")
    required = {"data", "indices", "indptr"}
    if not required <= set(group.keys()):
        raise TypeError(f"matrix {matrix_name!r} lacks {sorted(required)}")
    encoding = group.attrs.get("encoding-type", "")
    if isinstance(encoding, bytes):
        encoding = encoding.decode()
    if encoding not in {"", "csr_matrix"}:
        raise TypeError(f"only cell-major CSR h5ad matrices are supported, found {encoding!r}")
    return group


def sparse_shape(group: h5py.Group) -> tuple[int, int]:
    return tuple(int(x) for x in group.attrs["shape"])


def read_csr_rows(group: h5py.Group, start: int, stop: int, n_cols: int) -> sp.csr_matrix:
    pointers = group["indptr"][start : stop + 1]
    offset, end = int(pointers[0]), int(pointers[-1])
    return sp.csr_matrix(
        (
            group["data"][offset:end].astype(np.float64, copy=False),
            group["indices"][offset:end],
            pointers - offset,
        ),
        shape=(stop - start, n_cols),
    )


def h5_gene_names(h5: h5py.File, matrix_name: str) -> np.ndarray:
    var = h5["raw/var"] if matrix_name == "raw.X" else h5["var"]
    for column in ("feature_name", "gene_name", "gene_symbols", "symbol"):
        if column in var:
            return read_h5_string(var[column])
    key = var.attrs.get("_index", "_index")
    if isinstance(key, bytes):
        key = key.decode()
    return read_h5_string(var[key])


def transform_ids(values: pd.Series | np.ndarray, transform: str | None) -> np.ndarray:
    result = pd.Series(values, dtype="string").fillna("").astype(str)
    if transform in (None, "identity"):
        pass
    elif transform == "hyphen_to_underscore":
        result = result.str.replace("-", "_", regex=False)
    elif transform == "sample_prefix_hyphen_to_underscore":
        result = "sample_" + result.str.replace("-", "_", regex=False)
    else:
        raise ValueError(f"unknown sample transform {transform!r}")
    return result.to_numpy(dtype=str)


def load_manifest(path: str | Path) -> pd.DataFrame:
    frame = pd.read_csv(repo_path(path), dtype={"sample": str})
    if "sample" not in frame:
        raise ValueError(f"cohort manifest {path} has no sample column")
    if frame["sample"].duplicated().any():
        raise ValueError(f"cohort manifest {path} has duplicate samples")
    return frame


def map_labels(values: np.ndarray, mapping: dict[str, str] | None) -> np.ndarray:
    values = np.asarray(values).astype(str)
    if not mapping:
        return values
    return pd.Series(values).map(mapping).fillna("").to_numpy(dtype=str)


def matrix_market_header(path: str | Path) -> tuple[int, int, int]:
    with gzip.open(repo_path(path), "rt") as handle:
        line = handle.readline()
        while line.startswith("%"):
            line = handle.readline()
    return tuple(int(x) for x in line.split())


def iter_matrix_market(path: str | Path, chunk_nnz: int) -> Iterator[tuple[np.ndarray, np.ndarray, np.ndarray]]:
    with gzip.open(repo_path(path), "rt") as handle:
        line = handle.readline()
        while line.startswith("%"):
            line = handle.readline()
        reader = pd.read_csv(
            handle,
            sep=r"\s+",
            header=None,
            names=["row", "col", "value"],
            dtype={"row": np.int32, "col": np.int32, "value": np.float64},
            chunksize=chunk_nnz,
        )
        for chunk in reader:
            yield (
                chunk["row"].to_numpy() - 1,
                chunk["col"].to_numpy() - 1,
                chunk["value"].to_numpy(),
            )
