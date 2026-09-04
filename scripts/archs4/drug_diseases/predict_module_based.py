#!/usr/bin/env python3
# Module-based (LV-space) drug-disease prediction, one compendium at a time.
# Scores drug-disease pairs as -drug^T . disease in the given compendium's CLAMP
# LV space, at several top-N-LVs thresholds. Writes HDF5 prediction files under
# {output-dir}/lincs/predictions/dotprod_neg/ via drug_disease_utils.predict_dotprod_neg.

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from drug_disease_utils import predict_dotprod_neg


def parse_thresholds(raw: str) -> list[int | None]:
    return [None if tok == "all" else int(tok) for tok in raw.split(",")]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--compendium", required=True)
    parser.add_argument("--lincs-projection", required=True, type=Path)
    parser.add_argument("--spredixcan-proj-dir", required=True, type=Path)
    parser.add_argument("--gold-standard", required=True, type=Path)
    parser.add_argument("--ukb-efo", required=True, type=Path)
    parser.add_argument("--efo-xrefs", required=True, type=Path)
    parser.add_argument("--do-xrefs", required=True, type=Path)
    parser.add_argument("--n-top-lvs", required=True, help="comma-separated ints, or 'all' for no thresholding")
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()

    prediction_method = f"module_based_{args.compendium}"
    n_top_lvs_list = parse_thresholds(args.n_top_lvs)

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

    lincs_projection = pd.read_pickle(args.lincs_projection)
    print(f"LINCS projection shape: {lincs_projection.shape}")
    assert not lincs_projection.isna().any().any()

    spredixcan_file_list = sorted(
        args.spredixcan_proj_dir.glob(f"spredixcan-*-projection-{args.compendium}.pkl")
    )
    n_files = len(spredixcan_file_list)
    assert n_files == 49, f"Expected 49 tissue files, found {n_files}"

    for spredixcan_file in spredixcan_file_list:
        print(spredixcan_file.name)
        tissue_proj = pd.read_pickle(spredixcan_file)
        print(f"  shape: {tissue_proj.shape}")
        assert tissue_proj.index.equals(lincs_projection.index)

        for ntc in n_top_lvs_list:
            predict_dotprod_neg(
                lincs_projection,
                spredixcan_file,
                tissue_proj,
                output_predictions_dir,
                prediction_method,
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
