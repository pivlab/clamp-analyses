#!/usr/bin/env python3
# GTEx timming and RAM for all models, 3 seeds per each

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

# CoGAPS is intentionally absent: rule cogaps_gtex was run on GTEx and killed at
# its 7-day budget without converging.  See excluded_methods in
# workflow/config/runtime_benchmark_gtex.yaml.
METHODS = (
    "CLAMPfull",
    "CLAMPbase",
    "PLIER",
    "Flashier",
    "NMF",
    "PCA",
    "ICA",
    "GSSig",
    "MOFA-FLEX",
)

CHUNK = 64 << 20


def matrix_shape(path: Path) -> tuple[int, int]:
    with path.open("rb") as handle:
        header = handle.readline()
        n_samples = header.count(b",")
        n_genes = 0
        tail = b"\n"
        while chunk := handle.read(CHUNK):
            n_genes += chunk.count(b"\n")
            tail = chunk[-1:]
        if tail != b"\n":
            n_genes += 1
    return n_genes, n_samples


def read_shape(path: Path) -> tuple[int, int]:
    with path.open(newline="") as handle:
        row = next(csv.DictReader(handle))
    return int(row["n_genes"]), int(row["n_samples"])


def command_for(args: argparse.Namespace) -> tuple[list[str], list[str]]:
    scripts = Path("scripts/gtex")
    seed = str(args.seed)

    if args.method in ("CLAMPbase", "CLAMPfull"):
        return [
            "Rscript", str(scripts / "clamp.R"),
            "--fbm-filt", args.fbm_filt, "--svd-res", args.svd_res,
            "--k", args.k_rds, "--genes", args.genes, "--samples", args.samples,
            "--gmt", args.gmt, "--out-dir", args.output_dir,
            "--max-iter", str(args.clamp_max_iter), "--seed", seed,
            "--model", "base" if args.method == "CLAMPbase" else "full",
        ], [args.fbm_filt, args.svd_res]

    if args.method == "PLIER":
        return [
            "Rscript", str(scripts / "plier.R"),
            "--fbm-filt", args.fbm_filt, "--svd-res", args.svd_res,
            "--k", args.k_rds, "--genes", args.genes, "--samples", args.samples,
            "--gmt", args.gmt, "--out-dir", args.output_dir, "--seed", seed,
        ], [args.fbm_filt, args.svd_res]

    if args.method == "Flashier":
        return [
            "Rscript", str(scripts / "flashier.R"),
            "--df-gtex-fbm-filt", args.df_rds, "--k", args.k_rds,
            "--out-dir", args.output_dir,
            "--backfit-maxiter", str(args.flashier_backfit), "--seed", seed,
        ], [args.df_rds]

    if args.method == "GSSig":
        return [
            "Rscript", str(scripts / "gss.R"),
            "--dataset", args.dataset,
            "--df-gtex-fbm-filt", args.df_rds, "--k", args.k_rds,
            "--out-dir", args.output_dir,
            "--d-cluster", str(args.gss_d_cluster), "--seed", seed,
        ], [args.df_rds]

    if args.method in ("PCA", "NMF", "ICA"):
        return [
            sys.executable, str(scripts / "pca_nmf_ica.py"),
            "--df-gtex-fbm-filt", args.df_csv, "--k", args.k_csv,
            "--out-dir", args.output_dir, "--seed", seed,
            "--method", args.method, "--flat-output",
        ], [args.df_csv]

    if args.method == "MOFA-FLEX":
        return [
            sys.executable, str(scripts / "mofa_flex_prior.py"),
            "--df-gtex-fbm-filt", args.df_csv, "--k", args.k_csv,
            "--out-dir", args.output_dir,
            "--gmt", args.gmt,
            "--min-fraction", str(args.min_fraction),
            "--min-count", str(args.min_count),
            "--max-count", str(args.max_count),
            "--similarity-threshold", str(args.similarity_threshold),
            "--max-epochs", str(args.mofa_max_epochs), "--seed", seed,
            "--threads", str(args.threads),
        ], [args.df_csv]

    raise ValueError(f"Unsupported method: {args.method}")


def warm_page_cache(paths: list[str]) -> None:
    for path in paths:
        with open(path, "rb") as handle:
            while handle.read(CHUNK):
                pass


def directory_bytes(root: Path) -> int:
    return sum(p.stat().st_size for p in root.rglob("*") if p.is_file())


def write_timing(path: Path, row: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser()
    parser.add_argument("--probe-shape", action="store_true",
                        help="Write the matrix shape of --matrix-csv to --shape-csv "
                             "and exit.  Run once per benchmark, not once per fit.")
    parser.add_argument("--matrix-csv")
    parser.add_argument("--shape-csv")

    parser.add_argument("--dataset", default="GTEx")
    parser.add_argument("--method", choices=METHODS)
    parser.add_argument("--seed", type=int)
    parser.add_argument("--threads", type=int)

    parser.add_argument("--fbm-filt")
    parser.add_argument("--svd-res")
    parser.add_argument("--genes")
    parser.add_argument("--samples")
    parser.add_argument("--df-rds")
    parser.add_argument("--k-rds")
    parser.add_argument("--df-csv")
    parser.add_argument("--k-csv")

    parser.add_argument("--gmt")
    parser.add_argument("--output-dir")
    parser.add_argument("--timing-csv")
    parser.add_argument("--log-file")

    parser.add_argument("--clamp-max-iter", type=int, default=500)
    parser.add_argument("--flashier-backfit", type=int, default=20)
    parser.add_argument("--gss-d-cluster", type=int, default=4)
    parser.add_argument("--mofa-max-epochs", type=int, default=1000)
    parser.add_argument("--min-fraction", type=float, default=0.4)
    parser.add_argument("--min-count", type=int, default=40)
    parser.add_argument("--max-count", type=int, default=200)
    parser.add_argument("--similarity-threshold", type=float, default=0.8)
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.probe_shape:
        if not args.matrix_csv or not args.shape_csv:
            parser.error("--probe-shape requires --matrix-csv and --shape-csv")
        n_genes, n_samples = matrix_shape(Path(args.matrix_csv))
        out = Path(args.shape_csv)
        out.parent.mkdir(parents=True, exist_ok=True)
        with out.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["matrix", "n_genes", "n_samples"])
            writer.writerow([args.matrix_csv, n_genes, n_samples])
        print(f"[gtex] {args.matrix_csv}: {n_genes} genes x {n_samples} samples")
        return

    for flag in ("method", "seed", "threads", "output_dir", "timing_csv",
                 "log_file", "shape_csv"):
        if getattr(args, flag) is None:
            parser.error(f"the following arguments are required: --{flag.replace('_', '-')}")

    output_dir = Path(args.output_dir)
    timing_csv = Path(args.timing_csv)
    log_file = Path(args.log_file)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_file.parent.mkdir(parents=True, exist_ok=True)

    command, cache_inputs = command_for(args)

    env = os.environ.copy()
    for variable in (
        "OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
        "RCPP_PARALLEL_NUM_THREADS",
    ):
        env[variable] = str(args.threads)
    env["CUDA_VISIBLE_DEVICES"] = ""

    warm_page_cache(cache_inputs)

    started_at = datetime.now(timezone.utc)
    started = time.perf_counter()
    with log_file.open("w") as log:
        log.write(f"command: {shlex.join(command)}\n")
        log.write(f"threads: {args.threads}\n")
        log.flush()
        process = subprocess.Popen(
            command, cwd=Path.cwd(), env=env, stdout=log, stderr=subprocess.STDOUT
        )
        _pid, status, usage = os.wait4(process.pid, 0)
    elapsed_seconds = time.perf_counter() - started
    finished_at = datetime.now(timezone.utc)

    returncode = os.waitstatus_to_exitcode(status)
    process.returncode = returncode
    n_genes, n_samples = read_shape(Path(args.shape_csv))
    term_signal = -returncode if returncode < 0 else 0
    if returncode == 0:
        run_status = "success"
    elif term_signal:
        run_status = "killed"
    else:
        run_status = "failed"

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
        "status": run_status,
        "exit_code": returncode,
        "peak_rss_mb": f"{usage.ru_maxrss / 1024:.3f}",
        "peak_rss_source": "os.wait4:ru_maxrss",
        "user_cpu_seconds": f"{usage.ru_utime:.3f}",
        "system_cpu_seconds": f"{usage.ru_stime:.3f}",
        "cpu_utilization": f"{(usage.ru_utime + usage.ru_stime) / elapsed_seconds:.4f}",
        "started_epoch": f"{started_at.timestamp():.6f}",
        "finished_epoch": f"{finished_at.timestamp():.6f}",
        "term_signal": term_signal,
        "output_bytes": directory_bytes(output_dir),
    }
    write_timing(timing_csv, row)
    if returncode != 0:
        raise SystemExit(returncode)


if __name__ == "__main__":
    main()
