#!/usr/bin/env python3
"""Drop UK Biobank phenotypes in excluded Category-tree branches from a
combined phenoplier GLS trait table.

Runs inside the phenoplier-cli-neo conda env: category lookup uses
phenoplier's own `phenoplier.entity.Trait`, which derives a phenotype's
category from the UK Biobank data dictionary shipped in the linked reference
data bundle -- not reimplemented here. Administrative/behavioral UKB
categories (employment, food preferences, household, education, transport,
leisure, questionnaire admin fields, ...) are not health phenotypes and
otherwise inflate or deflate cross-model "traits recovered" comparisons.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
from pathlib import Path

import pandas as pd


@lru_cache(maxsize=None)
def phenotype_category(full_code: str) -> str | None:
    from phenoplier.entity import Trait

    try:
        return Trait.get_trait(full_code=full_code).category
    except Exception:
        return None


def is_excluded(category: str | None, excluded_terms: list[str]) -> bool:
    if category is None:
        return False
    category_lower = category.lower()
    return any(term in category_lower for term in excluded_terms)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--excluded-log", required=True, type=Path)
    parser.add_argument("--exclude", nargs="+", required=True)
    parser.add_argument("--phenotype-column", default="phenotype")
    args = parser.parse_args()

    excluded_terms = [term.lower() for term in args.exclude]

    df = pd.read_csv(args.input)
    df["ukb_category"] = df[args.phenotype_column].map(phenotype_category)
    excluded_mask = df["ukb_category"].apply(lambda c: is_excluded(c, excluded_terms))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.loc[~excluded_mask].to_csv(args.output, index=False)

    args.excluded_log.parent.mkdir(parents=True, exist_ok=True)
    (
        df.loc[excluded_mask, [args.phenotype_column, "ukb_category"]]
        .drop_duplicates()
        .to_csv(args.excluded_log, index=False)
    )


if __name__ == "__main__":
    main()
