#!/usr/bin/env python3
"""Concatenate every per-model phenoplier GLS summary into one long table."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

SEGMENT_RE = re.compile(r"^(rs|k|seed)(\d+)$")


def parse_model_key(model_key: str) -> dict:
    parts = model_key.split("/")
    fields = {"group": parts[0], "dataset": None, "fraction": None, "k": None, "seed": None}
    for part in parts[1:-1]:
        match = SEGMENT_RE.match(part)
        if not match:
            fields["dataset"] = part
            continue
        kind, value = match.groups()
        fields[{"rs": "fraction", "k": "k", "seed": "seed"}[kind]] = int(value)
    fields["dataset"] = fields["dataset"] or fields["group"]
    fields["model_type"] = parts[-1]
    return fields


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    frames = []
    for summary_path in sorted(args.root.glob("*/**/gls-summary-phenomexcan.tsv.gz")):
        model_key = str(summary_path.parent.relative_to(args.root))
        df = pd.read_csv(summary_path, sep="\t")
        for field, value in parse_model_key(model_key).items():
            df[field] = value
        df["model_key"] = model_key
        frames.append(df)

    if not frames:
        raise SystemExit(f"No gls-summary-phenomexcan.tsv.gz files found under {args.root}")

    combined = pd.concat(frames, ignore_index=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()
