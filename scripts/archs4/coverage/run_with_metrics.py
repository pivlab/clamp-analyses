#!/usr/bin/env python3

# Runs a model command and records its resource metrics

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import resource
import signal
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path


# Returns the current UTC timestamp.
def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


# Calculates a file checksum.
def file_hash(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


# Looks up a Slurm job state when available.
def slurm_state(job_id: str | None) -> str | None:
    if not job_id:
        return None
    try:
        result = subprocess.run(
            ["scontrol", "show", "job", "-o", job_id],
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    match = re.search(r"\bJobState=(\S+)", result.stdout)
    return match.group(1) if match else None


# Writes JSON without leaving a partial file.
def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(descriptor, "w") as handle:
            json.dump(payload, handle, indent=2)
            handle.write("\n")
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


# Runs the requested command and stores its metrics.
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--metrics-dir", required=True, type=Path)
    parser.add_argument("--requested-mem-mb", required=True, type=int)
    parser.add_argument("--initial-mem-mb", required=True, type=int)
    parser.add_argument("--requested-cpus", required=True, type=int)
    parser.add_argument("--requested-runtime-min", required=True, type=int)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--prior", type=Path)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        parser.error("a command is required after --")

    ratio = args.requested_mem_mb / args.initial_mem_mb
    attempt = 1 if ratio <= 1.01 else 2 if ratio <= 1.51 else 3
    output = args.metrics_dir / f"attempt{attempt}_mem{args.requested_mem_mb}.json"
    slurm_job_id = os.environ.get("SLURM_JOB_ID")
    started = utc_now()
    started_monotonic = time.monotonic()
    child = subprocess.Popen(command)

    def forward_signal(signum: int, _frame: object) -> None:
        if child.poll() is None:
            child.send_signal(signum)

    signal.signal(signal.SIGTERM, forward_signal)
    signal.signal(signal.SIGINT, forward_signal)
    exit_code = child.wait()
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    payload = {
        "schema_version": 1,
        "attempt": attempt,
        "status": "complete" if exit_code == 0 else "failed",
        "exit_code": exit_code,
        "started_utc": started,
        "finished_utc": utc_now(),
        "elapsed_seconds": time.monotonic() - started_monotonic,
        "max_rss_kb": usage.ru_maxrss,
        "requested": {
            "cpus": args.requested_cpus,
            "mem_mb": args.requested_mem_mb,
            "runtime_min": args.requested_runtime_min,
        },
        "rng_seed": args.seed,
        "prior": None if args.prior is None else {
            "path": str(args.prior.resolve()),
            "md5": file_hash(args.prior, "md5"),
            "sha256": file_hash(args.prior, "sha256"),
        },
        "slurm_job_id": slurm_job_id,
        "slurm_state_at_recording": slurm_state(slurm_job_id),
        "command": command,
    }
    atomic_json(output, payload)
    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()
