#!/usr/bin/env python3
"""Run one pseudobulk factorization and record auditable wall-clock timing."""

from __future__ import annotations

import argparse
import csv
import os
import shlex
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path


METHODS = (
    "CLAMPfull",
    "CLAMPbase",
    "PLIER",
    "CoGAPS",
    "Flashier",
    "NMF",
    "PCA",
    "ICA",
    "GSSig",
    "MOFA-FLEX",
)


def matrix_shape(path: Path) -> tuple[int, int]:
    """Return genes x samples for a CSV matrix with row names in column one."""
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        n_samples = len(header) - 1
        n_genes = sum(1 for _ in reader)
    return n_genes, n_samples


def command_for(args: argparse.Namespace) -> tuple[list[str], Path]:
    scripts = Path("scripts/pseudobulk")
    out = Path(args.output_dir)
    seed = str(args.seed)

    if args.method in ("CLAMPbase", "CLAMPfull"):
        command = [
            "Rscript", str(scripts / "clamp.R"),
            "--norm", args.norm,
            "--k", args.k,
            "--gmt", args.gmt,
            "--base-dir", args.output_dir,
            "--full-dir", args.output_dir,
            "--max-iter", str(args.clamp_max_iter),
            "--seed", seed,
            "--model", "base" if args.method == "CLAMPbase" else "full",
        ]
        return command, Path(args.norm)

    if args.method == "PLIER":
        return [
            "Rscript", str(scripts / "plier.R"),
            "--norm", args.norm, "--k", args.k, "--gmt", args.gmt,
            "--out-dir", args.output_dir, "--seed", seed,
        ], Path(args.norm)

    if args.method in ("PCA", "NMF", "ICA"):
        return [
            sys.executable, str(scripts / "pca_nmf_ica.py"),
            "--norm", args.norm, "--k", args.k,
            "--models-dir", args.output_dir, "--seed", seed,
            "--method", args.method, "--flat-output",
        ], Path(args.norm)

    if args.method == "Flashier":
        return [
            "Rscript", str(scripts / "flashier.R"),
            "--norm", args.norm, "--k", args.k,
            "--out-dir", args.output_dir,
            "--backfit-maxiter", str(args.flashier_backfit), "--seed", seed,
        ], Path(args.norm)

    if args.method == "MOFA-FLEX":
        return [
            sys.executable, str(scripts / "mofa_flex.py"),
            "--norm", args.norm, "--k", args.k, "--gmt", args.gmt,
            "--out-dir", args.output_dir,
            "--max-epochs", str(args.mofa_max_epochs), "--seed", seed,
            "--threads", str(args.threads),
        ], Path(args.norm)

    if args.method == "GSSig":
        return [
            "Rscript", str(scripts / "gss.R"),
            "--dataset", args.dataset, "--norm", args.norm, "--k", args.k,
            "--out-dir", args.output_dir, "--seed", seed,
        ], Path(args.norm)

    if args.method == "CoGAPS":
        return [
            "Rscript", str(scripts / "cogaps.R"),
            "--cpm-filt", args.cpm_filt, "--k", args.k,
            "--out-dir", args.output_dir,
            "--n-iterations", str(args.cogaps_iterations),
            "--n-threads", str(args.threads),
            "--parallel-workers", str(args.threads), "--seed", seed,
        ], Path(args.cpm_filt)

    raise ValueError(f"Unsupported method: {args.method}")


def write_timing(path: Path, row: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--method", required=True, choices=METHODS)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--threads", required=True, type=int)
    parser.add_argument("--norm", required=True)
    parser.add_argument("--cpm-filt", required=True)
    parser.add_argument("--k", required=True)
    parser.add_argument("--gmt", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--timing-csv", required=True)
    parser.add_argument("--log-file", required=True)
    parser.add_argument("--clamp-max-iter", type=int, default=500)
    parser.add_argument("--flashier-backfit", type=int, default=10)
    parser.add_argument("--mofa-max-epochs", type=int, default=200)
    parser.add_argument("--cogaps-iterations", type=int, default=5000)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    timing_csv = Path(args.timing_csv)
    log_file = Path(args.log_file)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_file.parent.mkdir(parents=True, exist_ok=True)

    command, matrix_path = command_for(args)
    env = os.environ.copy()
    for variable in (
        "OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
        "RCPP_PARALLEL_NUM_THREADS",
    ):
        env[variable] = str(args.threads)
    env["CUDA_VISIBLE_DEVICES"] = ""

    started_at = datetime.now(timezone.utc)
    started = time.perf_counter()
    with log_file.open("w") as log:
        log.write(f"command: {shlex.join(command)}\n")
        log.write(f"threads: {args.threads}\n")
        log.flush()
        completed = subprocess.run(
            command,
            cwd=Path.cwd(),
            env=env,
            stdout=log,
            stderr=subprocess.STDOUT,
            check=False,
        )
    elapsed_seconds = time.perf_counter() - started
    finished_at = datetime.now(timezone.utc)
    n_genes, n_samples = matrix_shape(matrix_path)

    row = {
        "dataset": args.dataset,
        "method": args.method,
        "seed": args.seed,
        "threads": args.threads,
        "n_genes": n_genes,
        "n_samples": n_samples,
        "elapsed_seconds": f"{elapsed_seconds:.9f}",
        "started_at": started_at.isoformat(),
        "finished_at": finished_at.isoformat(),
        "command": shlex.join(command),
        "status": "success" if completed.returncode == 0 else "failed",
        "exit_code": completed.returncode,
    }
    write_timing(timing_csv, row)
    if completed.returncode != 0:
        raise SystemExit(completed.returncode)


if __name__ == "__main__":
    main()
