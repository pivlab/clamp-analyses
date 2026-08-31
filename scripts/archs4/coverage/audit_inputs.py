#!/usr/bin/env python3
# Audits coverage cells, models, ORA results, and published BP models.
# Reports missing or inconsistent artifacts without changing them.

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import tempfile
from pathlib import Path

import yaml

CELL_ARTIFACTS = ("subsample_info.rds", "svd.rds", "CLAMP_K.rds", "CLAMPbase.rds")
MODEL_ARTIFACTS = ("CLAMPfull_bp.rds", "B.csv", "Z.csv", "summary.csv", "manifest.json")
ORA_ARTIFACTS = ("summary.csv", "enrichment.csv.gz")


# Reads RDS cell provenance through R.
# Reads coverage-cell provenance from RDS and HDF5 files.
def cell_provenance(args: argparse.Namespace, root: Path) -> dict[Path, dict]:
    r_code = r'''
suppressPackageStartupMessages(library(hdf5r))
args <- commandArgs(trailingOnly = TRUE)
h5_path <- args[[1L]]
input_root <- args[[2L]]
output <- args[[3L]]
sc_max <- as.numeric(args[[4L]])
`%||%` <- function(x, y) if (is.null(x)) y else x
h5 <- H5File$new(h5_path, mode = "r")
raw_samples <- h5[["/meta/samples/geo_accession"]]$read()
sc_probability <- h5[["/meta/samples/singlecellprobability"]]$read()
h5$close_all()
filtered_samples <- raw_samples[which(sc_probability < sc_max)]
paths <- sort(Sys.glob(file.path(input_root, "*", "study_coverage*", "subsample_info.rds")))
rows <- lapply(paths, function(info_path) {
  info <- readRDS(info_path)
  sample_idx <- as.integer(info$sample_idx)
  sample_names <- as.character(info$sample_names)
  idx_in_filtered <- length(sample_idx) == length(sample_names) && all(!is.na(sample_idx) & sample_idx >= 1L & sample_idx <= length(filtered_samples))
  idx_in_raw <- length(sample_idx) == length(sample_names) && all(!is.na(sample_idx) & sample_idx >= 1L & sample_idx <= length(raw_samples))
  status <- if (idx_in_filtered && identical(filtered_samples[sample_idx], sample_names)) "filtered_aligned" else if (idx_in_raw && identical(raw_samples[sample_idx], sample_names)) "raw_indexed_misaligned" else "inconsistent"
  data.frame(path = normalizePath(dirname(info_path)), sampling_provenance = status, declared_samples = as.integer(info$n_samples %||% length(sample_names)), stringsAsFactors = FALSE)
})
out <- if (length(rows)) do.call(rbind, rows) else data.frame(path = character(), sampling_provenance = character(), declared_samples = integer())
write.csv(out, output, row.names = FALSE)
'''
    with tempfile.TemporaryDirectory(prefix="coverage-audit-") as tmp:
        out = Path(tmp) / "cell_provenance.csv"
        subprocess.run(
            [
                "Rscript", "-e", r_code,
                str(root / args.raw_h5), str(root / args.input_root), str(out),
                str(args.single_cell_probability_max),
            ],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        with out.open(newline="") as handle:
            return {Path(row["path"]).resolve(): row for row in csv.DictReader(handle)}


# Audits coverage-cell artifacts and sample provenance.
def rows_for_cells(cfg: dict, root: Path, provenance: dict[Path, dict]) -> list[dict]:
    cov = cfg["archs4"]["coverage"]
    input_root = Path(cov["input_root"])
    rows = []
    for fraction, level_dir in cov["levels"].items():
        for seed in cov["seeds"]:
            cell = root / input_root / level_dir / f"study_coverage_rs{fraction}_seed_{seed}"
            wanted = list(CELL_ARTIFACTS)
            if int(fraction) != 100:
                wanted.append("fbm_subsampled.bk")
            missing = [name for name in wanted if not (cell / name).exists()]
            audit = provenance.get(cell.resolve(), {})
            sampling = audit.get("sampling_provenance", "not_audited")
            status = "ok" if cell.is_dir() and not missing and sampling == "filtered_aligned" else "incomplete"
            if cell.is_dir() and not missing and sampling == "raw_indexed_misaligned":
                status = "provenance_mismatch"
            rows.append({
                "kind": "cell",
                "dataset": "archs4",
                "fraction": fraction,
                "seed": seed,
                "path": str(cell),
                "exists": cell.is_dir(),
                "missing": ";".join(missing),
                "sampling_provenance": sampling,
                "declared_samples": audit.get("declared_samples", ""),
                "status": status,
            })
    return rows


# Audits fitted GO:BP model artifacts.
def rows_for_models(cfg: dict, root: Path) -> list[dict]:
    cov = cfg["archs4"]["coverage"]
    model_root = Path(cov["model_root"])
    model_name = cov["model_name"]
    datasets = ["archs4"] + list(cov["comparators"])
    rows = []
    for dataset in datasets:
        fractions = list(cov["levels"]) if dataset == "archs4" else ["100"]
        for fraction in fractions:
            for seed in cov["seeds"]:
                base = root / model_root / dataset / f"rs{fraction}" / f"seed{seed}"
                model_dir = base / model_name
                missing = [n for n in MODEL_ARTIFACTS if not (model_dir / n).is_file()]
                validated = base / "validated.json"
                stale = ""
                if validated.is_file():
                    try:
                        recorded = json.loads(validated.read_text()).get("model_dir", "")
                    except json.JSONDecodeError:
                        recorded = ""
                    if recorded and Path(recorded).resolve() != model_dir.resolve():
                        stale = recorded
                rows.append({
                    "kind": "model",
                    "dataset": dataset,
                    "fraction": fraction,
                    "seed": seed,
                    "path": str(model_dir),
                    "exists": model_dir.is_dir(),
                    "missing": ";".join(missing),
                    "validated": validated.is_file(),
                    "stale_validated_path": stale,
                    "status": "ok" if model_dir.is_dir() and not missing else "incomplete",
                })
    return rows


# Audits completed ORA result directories.
def rows_for_ora(cfg: dict, root: Path) -> list[dict]:
    cov = cfg["archs4"]["coverage"]
    ora_root = root / cov["ora_root"]
    databases = list(cov["databases"])
    rows = []
    if not ora_root.is_dir():
        return rows
    for summary in sorted(ora_root.rglob("summary.csv")):
        rel = summary.parent.relative_to(ora_root).parts
        if len(rel) != 5:
            continue
        dataset, fraction, seed, model, database = rel
        missing = [n for n in ORA_ARTIFACTS if not (summary.parent / n).is_file()]
        rows.append({
            "kind": "ora",
            "dataset": dataset,
            "fraction": fraction.removeprefix("rs"),
            "seed": seed.removeprefix("seed"),
            "model": model,
            "database": database,
            "path": str(summary.parent),
            "exists": True,
            "missing": ";".join(missing),
            "known_database": database in databases,
            "status": "ok" if not missing else "incomplete",
        })
    return rows


# Audits the published seed-free models.
def rows_for_bp(cfg: dict, root: Path) -> list[dict]:
    cov = cfg["archs4"]["coverage"]
    bp_root = root / cov["bp_models_root"]
    rows = []
    for dataset in cov["bp_model_datasets"]:
        out = bp_root / dataset
        payload = out / "CLAMPfull_bp"
        wanted = ("CLAMPfull_bp.rds", "B.csv", "Z.csv", "summary.csv")
        missing = [n for n in wanted if not (payload / n).is_file()]
        if not (out / "manifest.json").is_file():
            missing.append("manifest.json")
        rows.append({
            "kind": "bp_model",
            "dataset": dataset,
            "path": str(out),
            "exists": out.is_dir(),
            "missing": ";".join(missing),
            "status": "ok" if out.is_dir() and not missing else "incomplete",
        })
    return rows


# Writes the combined audit report.
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--repo-root", type=Path, default=Path("."))
    parser.add_argument("--csv-out", required=True, type=Path)
    parser.add_argument("--json-out", required=True, type=Path)
    parser.add_argument(
        "--strict",
        action="store_true",
        help="exit non-zero when anything is missing (default: report and exit 0)",
    )
    args = parser.parse_args()

    cfg = yaml.safe_load(args.config.read_text())
    root = args.repo_root
    args.raw_h5 = cfg["archs4"]["raw_h5"]
    args.input_root = cfg["archs4"]["coverage"]["input_root"]
    args.single_cell_probability_max = cfg["archs4"]["preprocess"]["single_cell_probability_max"]
    provenance = cell_provenance(args, root)

    rows = (
        rows_for_cells(cfg, root, provenance)
        + rows_for_models(cfg, root)
        + rows_for_ora(cfg, root)
        + rows_for_bp(cfg, root)
    )

    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)

    args.csv_out.parent.mkdir(parents=True, exist_ok=True)
    with args.csv_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    summary: dict[str, dict[str, int]] = {}
    for row in rows:
        bucket = summary.setdefault(row["kind"], {})
        bucket[row["status"]] = bucket.get(row["status"], 0) + 1
    stale = [r["path"] for r in rows if r.get("stale_validated_path")]
    report = {
        "counts": summary,
        "incomplete": [r["path"] for r in rows if r["status"] != "ok"],
        "stale_validated_paths": stale,
    }
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(report, indent=2) + "\n")

    for kind, counts in summary.items():
        detail = " ".join(f"{status}={count:4d}" for status, count in sorted(counts.items()))
        print(f"{kind:9s} {detail}")
    if stale:
        print(f"\n{len(stale)} validated.json file(s) record a stale model_dir path")
    for path in report["incomplete"]:
        print(f"  incomplete: {path}")

    if args.strict and report["incomplete"]:
        raise SystemExit(f"{len(report['incomplete'])} incomplete entries")


if __name__ == "__main__":
    main()
