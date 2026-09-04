#!/usr/bin/env python3
# Single-gene-based drug-disease prediction baseline (no CLAMP projection). For
# each of the 49 raw (cleaned, but unprojected) S-PrediXcan tissue files, intersects
# genes with the raw LINCS matrix and scores drug-disease pairs as -drug^T . disease
# in gene space, at several top-N-genes thresholds. This step never touches a CLAMP
# model, so it is run once (not per compendium). Writes HDF5 prediction files under
# {output-dir}/lincs/predictions/dotprod_neg/ via drug_disease_utils.predict_dotprod_neg.

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from drug_disease_utils import predict_dotprod_neg

PREDICTION_METHOD = "gene_based"


def parse_thresholds(raw: str) -> list[int | None]:
    return [None if tok == "all" else int(tok) for tok in raw.split(",")]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spredixcan-raw-dir", required=True, type=Path)
    parser.add_argument("--lincs-file", required=True, type=Path)
    parser.add_argument("--gold-standard", required=True, type=Path)
    parser.add_argument("--ukb-efo", required=True, type=Path)
    parser.add_argument("--efo-xrefs", required=True, type=Path)
    parser.add_argument("--do-xrefs", required=True, type=Path)
    parser.add_argument("--n-top-genes", required=True, help="comma-separated ints, or 'all' for no thresholding")
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    n_top_genes_list = parse_thresholds(args.n_top_genes)

    output_predictions_dir = args.output_dir / "lincs" / "predictions"
    output_predictions_dir.mkdir(parents=True, exist_ok=True)

    gold_standard = pd.read_pickle(args.gold_standard)
    print(f"Gold standard shape: {gold_standard.shape}")
    doids_in_gold_standard = set(gold_standard["trait"])
    print(f"Unique DOIDs in gold standard: {len(doids_in_gold_standard)}")

    ukb_efo = pd.read_csv(args.ukb_efo, sep="\t", index_col="ukb_fullcode")
    ukb_efo.index = [idx.replace("-", "_", 1) for idx in ukb_efo.index]

    efo_xrefs = pd.read_csv(args.efo_xrefs, sep="\t")
    do_xrefs = pd.read_csv(args.do_xrefs, sep="\t")

    lincs_data = pd.read_pickle(args.lincs_file)
    print(f"LINCS shape: {lincs_data.shape}")

    spredixcan_file_list = sorted(
        f for f in args.spredixcan_raw_dir.glob("*.pkl") if f.name.startswith("spredixcan-")
    )
    n_files = len(spredixcan_file_list)
    assert n_files == 49, f"Expected 49 tissue files, found {n_files}"

    for spredixcan_file in spredixcan_file_list:
        print(spredixcan_file.name)
        tissue_data = pd.read_pickle(spredixcan_file)

        common_genes = tissue_data.index.intersection(lincs_data.index)
        tissue_data_common = tissue_data.loc[common_genes]
        lincs_common = lincs_data.loc[common_genes]
        print(f"  shape: {tissue_data_common.shape} ({len(common_genes)} common genes)")

        for ntc in n_top_genes_list:
            predict_dotprod_neg(
                lincs_common,
                spredixcan_file,
                tissue_data_common,
                output_predictions_dir,
                PREDICTION_METHOD,
                doids_in_gold_standard,
                ukb_efo,
                efo_xrefs,
                do_xrefs,
                n_top_conditions=ntc,
                use_abs=True,
            )
        print("")


if __name__ == "__main__":
    main()
