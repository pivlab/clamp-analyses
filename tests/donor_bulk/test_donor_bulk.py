"""Synthetic regression tests for donor-bulk aggregation and projection semantics."""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "donor_bulk"))
sys.path.insert(0, str(ROOT / "scripts" / "pseudobulk"))

from build_donor_bulk import (  # noqa: E402
    aggregate_h5,
    aggregate_mtx,
    gene_aggregator,
    select_cells,
)
from single_cell_projection import raw_feature_weights  # noqa: E402
from common import filter_analysis_cell_types  # noqa: E402


COUNTS = np.asarray(
    [
        [1, 2, 3, 0],
        [2, 0, 1, 1],
        [0, 3, 0, 2],
        [4, 0, 2, 0],
        [7, 1, 0, 1],  # excluded cell
        [1, 1, 1, 1],
    ],
    dtype=np.int64,
)
GENES = np.asarray(["G1", "G1", "G2", "G3"])
CELL_IDS = np.asarray([f"cell{i}" for i in range(len(COUNTS))])
DONORS = np.asarray(["D1", "D1", "D2", "D2", "D2", "D1"])
LABELS = np.asarray(["A", "B", "A", "B", "Drop", "A"])
MANIFEST = pd.DataFrame({"sample": ["D1", "D2"], "batch": ["x", "y"]})


def _write_sparse_group(parent: h5py.File, key: str, matrix: np.ndarray) -> None:
    csr = sp.csr_matrix(matrix)
    group = parent.create_group(key)
    group.attrs["encoding-type"] = "csr_matrix"
    group.attrs["shape"] = matrix.shape
    group.create_dataset("data", data=csr.data)
    group.create_dataset("indices", data=csr.indices)
    group.create_dataset("indptr", data=csr.indptr)


def _write_h5ad(path: Path) -> None:
    with h5py.File(path, "w") as h5:
        _write_sparse_group(h5, "X", COUNTS)
        obs = h5.create_group("obs")
        obs.attrs["_index"] = "_index"
        string = h5py.string_dtype("utf-8")
        obs.create_dataset("_index", data=CELL_IDS.astype(object), dtype=string)
        obs.create_dataset("donor", data=DONORS.astype(object), dtype=string)
        obs.create_dataset("label", data=LABELS.astype(object), dtype=string)
        var = h5.create_group("var")
        var.attrs["_index"] = "_index"
        var.create_dataset("_index", data=GENES.astype(object), dtype=string)


def _write_mtx_fixture(directory: Path) -> dict:
    matrix = directory / "matrix.mtx.gz"
    with gzip.open(matrix, "wt") as handle:
        handle.write("%%MatrixMarket matrix coordinate integer general\n")
        handle.write("% synthetic donor-bulk fixture\n")
        handle.write(f"{COUNTS.shape[1]} {COUNTS.shape[0]} {np.count_nonzero(COUNTS)}\n")
        for cell in range(COUNTS.shape[0]):
            for gene in range(COUNTS.shape[1]):
                if COUNTS[cell, gene]:
                    handle.write(f"{gene + 1} {cell + 1} {COUNTS[cell, gene]}\n")
    features = directory / "features.tsv.gz"
    barcodes = directory / "barcodes.tsv.gz"
    with gzip.open(features, "wt") as handle:
        handle.write("\n".join(GENES) + "\n")
    with gzip.open(barcodes, "wt") as handle:
        handle.write("\n".join(CELL_IDS) + "\n")
    metadata = directory / "meta.tsv"
    pd.DataFrame({"cell": CELL_IDS, "donor": DONORS, "label": LABELS}).to_csv(
        metadata, sep="\t", index=False
    )
    return {
        "kind": "matrix_market",
        "raw": str(matrix),
        "features": str(features),
        "barcodes": str(barcodes),
        "metadata": str(metadata),
        "cell_id_col": "cell",
        "sample_col": "donor",
        "sample_transform": "identity",
        "truth_v0_col": "label",
        "exclude": {"label": ["Drop"]},
    }


def _assert_result(result: tuple) -> None:
    sums, mean_cpm, genes, truth, donors, denominators, sampled = result[:7]
    keep = LABELS != "Drop"
    expected = np.vstack([COUNTS[keep & (DONORS == donor)].sum(axis=0) for donor in donors])
    np.testing.assert_array_equal(sums, expected)
    np.testing.assert_array_equal(denominators, [3, 2])
    np.testing.assert_allclose(truth.sum(axis=1), 1)
    np.testing.assert_array_equal(
        np.rint(truth.to_numpy() * denominators[:, None]).sum(axis=1), denominators
    )
    np.testing.assert_allclose(mean_cpm.sum(axis=1), 1e6)
    unique, aggregator = gene_aggregator(genes)
    assert list(unique) == ["G1", "G2", "G3"]
    np.testing.assert_array_equal(np.asarray(sums @ aggregator)[:, 0], expected[:, :2].sum(axis=1))
    assert sampled.shape[0] == 4


def test_h5ad_aggregation_exclusions_duplicates_and_denominators(tmp_path: Path) -> None:
    path = tmp_path / "fixture.h5ad"
    _write_h5ad(path)
    cfg = {
        "kind": "h5ad",
        "raw": str(path),
        "matrix": "X",
        "sample_col": "donor",
        "sample_transform": "identity",
        "truth_v0_col": "label",
        "exclude": {"label": ["Drop"]},
    }
    result = aggregate_h5("synthetic_h5", cfg, MANIFEST, 2, 4, 1, 123)
    _assert_result(result)


def test_matrix_market_matches_h5ad_policy(tmp_path: Path) -> None:
    cfg = _write_mtx_fixture(tmp_path)
    result = aggregate_mtx("synthetic_mtx", cfg, MANIFEST, 3, 4, 1, 123, tmp_path)
    _assert_result(result)


def test_sampling_is_deterministic_and_covers_each_type() -> None:
    eligible = LABELS != "Drop"
    first = select_cells(LABELS, eligible, target=2, minimum_per_type=1, seed=123)
    second = select_cells(LABELS, eligible, target=2, minimum_per_type=1, seed=123)
    np.testing.assert_array_equal(first, second)
    assert set(LABELS[first]) == {"A", "B"}
    assert not np.isin(first, np.flatnonzero(~eligible)).any()


def test_projection_duplicate_gene_policy_is_sum_consistent() -> None:
    model_genes = pd.Index(["G1", "G2"])
    model_weights = np.asarray([[2.0, 4.0], [3.0, 6.0]])
    raw_genes = np.asarray(["G1", "G1", "G2", "unused"])
    sum_weights, overlap = raw_feature_weights(
        raw_genes, model_genes, model_weights, duplicate_gene_policy="sum"
    )
    mean_weights, _ = raw_feature_weights(
        raw_genes, model_genes, model_weights, duplicate_gene_policy="mean"
    )
    assert overlap == 3
    np.testing.assert_array_equal(sum_weights[:2], np.vstack([model_weights[0]] * 2))
    np.testing.assert_array_equal(mean_weights[:2], np.vstack([model_weights[0] / 2] * 2))
    np.testing.assert_array_equal(sum_weights[2], model_weights[1])


def test_analysis_exclusions_are_explicit_and_do_not_renormalize() -> None:
    truth = pd.DataFrame(
        {"A": [0.60, 0.70], "Rare": [0.01, 0.02], "B": [0.39, 0.28]},
        index=["D1", "D2"],
    )
    config = {
        "cell_type_analysis": {
            "excluded_targets": {"synthetic": ["Rare"]}
        }
    }
    filtered = filter_analysis_cell_types(truth, config, "synthetic")
    assert list(filtered.columns) == ["A", "B"]
    np.testing.assert_allclose(filtered.sum(axis=1), [0.99, 0.98])
