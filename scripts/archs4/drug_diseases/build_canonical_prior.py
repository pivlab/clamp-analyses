#!/usr/bin/env python3
# Build the canonical CLAMP prior used by the drug-disease models.
# The generated GMT is deliberately self-contained and prefixes every term by
# its source collection. This prevents term-name collisions and makes the
# model manifest auditable without requiring an external database release.
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd


def read_gmt(path: Path, prefix: str):
    rows = []
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                rows.append((f"{prefix}_{fields[0]}", fields[1], list(dict.fromkeys(fields[2:]))))
    return rows


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--c2-all", required=True, type=Path)
    p.add_argument("--reactome", required=True, type=Path)
    p.add_argument("--cellmarker", required=True, type=Path)
    p.add_argument("--sheet", default="human")
    p.add_argument("--term-column", default="cell_name")
    p.add_argument("--gene-column", default="Symbol")
    p.add_argument("--output-gmt", required=True, type=Path)
    p.add_argument("--manifest", required=True, type=Path)
    args = p.parse_args()

    rows = read_gmt(args.c2_all, "C2ALL") + read_gmt(args.reactome, "REACTOME")
    markers = pd.read_excel(args.cellmarker, sheet_name=args.sheet)
    required = {args.term_column, args.gene_column}
    missing = required - set(markers.columns)
    if missing:
        raise ValueError(f"CellMarker columns missing: {sorted(missing)}")
    for term, group in markers.dropna(subset=list(required)).groupby(args.term_column, sort=True):
        genes = list(dict.fromkeys(group[args.gene_column].astype(str)))
        if genes:
            rows.append((f"CELLMARKER_{term}", "CellMarker_Human", genes))
    names = [row[0] for row in rows]
    if len(names) != len(set(names)):
        raise ValueError("Canonical prior contains duplicate prefixed set names")

    args.output_gmt.parent.mkdir(parents=True, exist_ok=True)
    with args.output_gmt.open("w") as handle:
        for name, desc, genes in rows:
            handle.write("\t".join([name, desc, *genes]) + "\n")
    digest = hashlib.md5(args.output_gmt.read_bytes()).hexdigest()
    manifest = {
        "schema_version": 1,
        "collections": {"C2ALL": sum(n.startswith("C2ALL_") for n in names),
                        "REACTOME": sum(n.startswith("REACTOME_") for n in names),
                        "CELLMARKER": sum(n.startswith("CELLMARKER_") for n in names)},
        "total_sets": len(rows), "md5": digest,
        "inputs": {k: str(getattr(args, k)) for k in ("c2_all", "reactome", "cellmarker")},
    }
    args.manifest.parent.mkdir(parents=True, exist_ok=True)
    args.manifest.write_text(json.dumps(manifest, indent=2) + "\n")
    if manifest["total_sets"] != 11107:
        raise ValueError(f"Expected 11107 canonical sets, found {manifest['total_sets']}")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
