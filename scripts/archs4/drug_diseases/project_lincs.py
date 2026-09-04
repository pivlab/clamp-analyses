#!/usr/bin/env python3
# Project the LINCS L1000 drug-perturbation signature matrix into a CLAMP LV space.
# Projects the drug x gene LINCS matrix through the given compendium's CLAMP model
# (--model-name; the canonical drug-diseases pipeline passes CLAMPfull_canonical)
# via CLAMP::projectCLAMP, producing a drug x LV matrix.

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from drug_disease_utils import prepare_clamp_projector, project_to_clamp


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-dir", required=True, type=Path)
    parser.add_argument("--model-name", required=True)
    parser.add_argument("--lincs-file", required=True, type=Path)
    parser.add_argument("--compendium", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    clamp_model_file = args.model_dir / f"{args.model_name}.rds"
    assert clamp_model_file.exists(), clamp_model_file

    lincs_data = pd.read_pickle(args.lincs_file)
    print(f"LINCS shape: {lincs_data.shape}")
    assert lincs_data.index.is_unique
    assert lincs_data.columns.is_unique
    assert not lincs_data.isna().any().any()

    CLAMP, clamp_sub, mapped_ensembl, lv_names = prepare_clamp_projector(clamp_model_file)

    lincs_projection = project_to_clamp(lincs_data, CLAMP, clamp_sub, mapped_ensembl, lv_names)
    print(f"LINCS projection shape: {lincs_projection.shape}")
    assert not lincs_projection.isna().any().any()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    lincs_projection.to_pickle(args.output)
    print(f"Saved to: {args.output}")


if __name__ == "__main__":
    main()
