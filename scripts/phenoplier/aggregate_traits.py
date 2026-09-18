#!/usr/bin/env python3
"""Concatenate per-model phenoplier GLS summaries into one long table.

Each input is a summary published by run_gls.sh, named the way
scripts/archs4/traits/aggregate_*_traits.R expect:

    cov_rs<fraction>_seed<seed>.tsv.gz          ARCHS4 coverage cell
    sat_rs<fraction>_k<k>_seed<seed>.tsv.gz     ARCHS4 saturation cell
    final_<dataset>.tsv.gz                      published full-data model
    canon_<dataset>.tsv.gz                      published canonical model

The model family (CLAMPfull_bp vs CLAMPbase) is not in the filename -- it is
the directory a summary lives in -- so the caller passes each family's files
under its own flag.

A summary is one row per LV x phenotype (~4 M rows for a 1,728-LV model), so
the table is streamed: each input is read in chunks and appended to the
gzipped output, never held in memory as a whole. Expect the output to be tens
of GB for the full model set.
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import pandas as pd

CHUNK_ROWS = 500_000

PATTERNS = {
    "coverage": re.compile(r"^cov_rs(?P<fraction>\d+)_seed(?P<seed>\d+)$"),
    "saturation": re.compile(r"^sat_rs(?P<fraction>\d+)_k(?P<k>\d+)_seed(?P<seed>\d+)$"),
    "final": re.compile(r"^final_(?P<dataset>[A-Za-z0-9]+)$"),
    "canonical": re.compile(r"^canon_(?P<dataset>[A-Za-z0-9]+)$"),
}


def parse_summary_name(path: Path) -> dict:
    stem = path.name[: -len(".tsv.gz")] if path.name.endswith(".tsv.gz") else path.stem
    for study, pattern in PATTERNS.items():
        match = pattern.match(stem)
        if match:
            fields = {"study": study, "dataset": "archs4", "fraction": None, "k": None, "seed": None}
            for key, value in match.groupdict().items():
                fields[key] = int(value) if key in {"fraction", "k", "seed"} else value
            return fields
    raise SystemExit(f"Unrecognised summary name: {path}")


def stream(paths: list[str], model: str, out, state: dict) -> None:
    for raw in paths:
        path = Path(raw)
        labels = parse_summary_name(path)
        labels["model"] = model
        labels["summary"] = str(path)
        rows = 0
        for chunk in pd.read_csv(path, sep="\t", low_memory=False, chunksize=CHUNK_ROWS):
            for field, value in labels.items():
                chunk[field] = value
            if state["columns"] is None:
                state["columns"] = list(chunk.columns)
            elif list(chunk.columns) != state["columns"]:
                raise SystemExit(f"Column mismatch in {path}: {list(chunk.columns)} != {state['columns']}")
            chunk.to_csv(out, index=False, header=not state["header_written"])
            state["header_written"] = True
            rows += len(chunk)
        state["files"] += 1
        state["rows"] += rows
        print(f"{path}: {rows} rows", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True, type=Path, help="output .csv.gz")
    parser.add_argument("--clampfull-bp", nargs="*", default=[], help="CLAMPfull_bp summaries")
    parser.add_argument("--clampbase", nargs="*", default=[], help="CLAMPbase summaries")
    parser.add_argument("--canonical", nargs="*", default=[], help="CLAMPfull_canonical summaries")
    args = parser.parse_args()

    if not (args.clampfull_bp or args.clampbase or args.canonical):
        raise SystemExit("No summaries given")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.out.with_name(args.out.name + ".tmp")
    state = {"columns": None, "header_written": False, "files": 0, "rows": 0}
    with gzip.open(tmp, "wt", newline="") as out:
        stream(args.clampfull_bp, "CLAMPfull_bp", out, state)
        stream(args.clampbase, "CLAMPbase", out, state)
        stream(args.canonical, "CLAMPfull_canonical", out, state)
    tmp.replace(args.out)
    print(f"{state['files']} summaries, {state['rows']} rows -> {args.out}")


if __name__ == "__main__":
    main()
