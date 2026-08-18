#!/usr/bin/env python3
"""Validate an atomically published coverage CLAMPfull model."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import subprocess
import tempfile
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def line_count(path: Path) -> int:
    with path.open("rb") as handle:
        return sum(1 for _ in handle)


def slurm_size_kb(value: str) -> int:
    match = re.fullmatch(r"([0-9.]+)([KMGT]?)", value)
    if not match:
        return 0
    factors = {"": 1 / 1024, "K": 1, "M": 1024, "G": 1024**2, "T": 1024**3}
    return int(float(match.group(1)) * factors[match.group(2)])


def slurm_accounting(job_id: str | None) -> dict | None:
    if not job_id or job_id == "NA":
        return None
    try:
        result = subprocess.run(
            [
                "sacct", "-j", job_id, "--noheader", "--parsable2",
                "--format=JobIDRaw,State,ExitCode,ElapsedRaw,MaxRSS,ReqMem,AllocCPUS",
            ],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    rows = [line.split("|") for line in result.stdout.splitlines() if line.strip()]
    if not rows:
        return None
    primary = rows[0]
    rss_values = [row[4] for row in rows if len(row) >= 5 and row[4]]
    max_rss = max(rss_values, key=slurm_size_kb, default="")
    return {
        "job_id": primary[0],
        "state": primary[1],
        "exit_code": primary[2],
        "elapsed_seconds": int(primary[3]) if primary[3].isdigit() else None,
        "max_rss": max_rss or None,
        "max_rss_kb": slurm_size_kb(max_rss) if max_rss else None,
        "requested_memory": primary[5] or None,
        "allocated_cpus": int(primary[6]) if primary[6].isdigit() else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-dir", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    model_dir = args.model_dir.resolve()
    manifest_path = model_dir / "manifest.json"
    required = {
        "completion marker": model_dir / "complete",
        "manifest": manifest_path,
        "model": model_dir / "CLAMPfull_bp.rds",
        "B": model_dir / "B.csv",
        "Z": model_dir / "Z.csv",
        "summary": model_dir / "summary.csv",
    }
    missing = [label for label, path in required.items() if not path.is_file()]
    if missing:
        raise SystemExit(f"Missing model artifacts: {', '.join(missing)}")
    empty = [label for label, path in required.items() if path.stat().st_size == 0]
    if empty:
        raise SystemExit(f"Empty model artifacts: {', '.join(empty)}")

    manifest = json.loads(manifest_path.read_text())
    if manifest.get("status") != "complete":
        raise SystemExit("Model manifest is not complete")
    genes = int(manifest["genes"])
    samples = int(manifest["samples"])
    latent_variables = int(manifest["latent_variables"])
    attempts_dir = model_dir.parent / "attempts"
    attempts = [
        json.loads(path.read_text())
        for path in sorted(attempts_dir.glob("attempt*.json"))
    ] if attempts_dir.is_dir() else []

    z_path = required["Z"]
    if line_count(z_path) != genes + 1:
        raise SystemExit("Z row count does not match the model manifest")
    with z_path.open(newline="") as handle:
        z_columns = len(next(csv.reader(handle)))
    if z_columns != latent_variables + 1:
        raise SystemExit("Z column count does not match the model manifest")

    b_path = required["B"]
    with b_path.open(newline="") as handle:
        b_columns = len(next(csv.reader(handle)))
    if b_columns != samples + 1:
        raise SystemExit("B column count does not match the model manifest")

    result = {
        "schema_version": 1,
        "status": "validated",
        "model_dir": str(model_dir),
        "dataset": manifest["dataset"],
        "fraction": manifest["fraction"],
        "seed_index": manifest["seed_index"],
        "genes": genes,
        "samples": samples,
        "latent_variables": latent_variables,
        "z_sha256": sha256(z_path),
        "manifest_sha256": sha256(manifest_path),
        "slurm_job_id": manifest.get("slurm_job_id"),
        "rng_seed": manifest.get("rng_seed"),
        "prior": manifest.get("prior"),
        "elapsed_seconds": manifest.get("elapsed_seconds"),
        "requested_resources": manifest.get("requested_resources"),
        "attempts": attempts,
        "slurm_accounting": slurm_accounting(manifest.get("slurm_job_id")),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(
        prefix=f".{args.output.name}.", dir=args.output.parent
    )
    try:
        with os.fdopen(descriptor, "w") as handle:
            json.dump(result, handle, indent=2)
            handle.write("\n")
        os.replace(temporary, args.output)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


if __name__ == "__main__":
    main()
