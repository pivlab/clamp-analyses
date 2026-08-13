#!/usr/bin/env python3
"""Route ready ARCHS4 ORA targets between a Slurm login host and this host.

The script contains no site-specific hostname or path.  Those are supplied at
runtime.  It is intentionally an outer dispatcher: each selected backend still
runs the declared Snakemake ORA target, while this process owns the assignment
manifest and prevents the same target from being launched twice.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import re
import shlex
import shutil
import subprocess
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


DATABASES = ("reactome", "canonical", "cellmarker")
VALIDATED_RE = re.compile(r"/archs4/rs(?P<fraction>\d+)/seed(?P<seed>\d+)/validated\.json$")
KEY_RE = re.compile(
    r"^archs4_rs(?P<fraction>\d+)_seed(?P<seed>\d+)_CLAMPfull_(?P<database>\w+)$"
)


@dataclass(frozen=True)
class Target:
    fraction: int
    seed: int
    database: str

    @classmethod
    def from_key(cls, key: str) -> "Target | None":
        match = KEY_RE.match(key)
        if not match or match["database"] not in DATABASES:
            return None
        return cls(int(match["fraction"]), int(match["seed"]), match["database"])

    @property
    def key(self) -> str:
        return f"archs4_rs{self.fraction}_seed{self.seed}_CLAMPfull_{self.database}"

    @property
    def model_leaf(self) -> str:
        return f"archs4/rs{self.fraction}/seed{self.seed}"

    @property
    def ora_leaf(self) -> str:
        return f"archs4/rs{self.fraction}/seed{self.seed}/CLAMPfull/{self.database}"


def run(
    command: list[str], *, capture: bool = False, cwd: str | None = None
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        command,
        check=True,
        text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.PIPE if capture else None,
        cwd=cwd,
    )


def ssh(remote: str, command: str, *, capture: bool = False) -> subprocess.CompletedProcess[str]:
    return run(["ssh", "-o", "BatchMode=yes", remote, command], capture=capture)


def remote_capacity(
    remote: str, min_cpus: int, min_mem_mb: int
) -> tuple[int, int, int, int, bool]:
    nodes = ssh(remote, "scontrol show node -o", capture=True).stdout.splitlines()
    free_cpus = 0
    free_mem = 0
    ora_capacity = 0
    for node in nodes:
        values = dict(re.findall(r"(CPUAlloc|CPUTot|RealMemory|AllocMem)=(\d+)", node))
        if values:
            node_cpus = int(values["CPUTot"]) - int(values["CPUAlloc"])
            node_mem = int(values["RealMemory"]) - int(values["AllocMem"])
            free_cpus += node_cpus
            free_mem += node_mem
            ora_capacity += min(node_cpus // min_cpus, node_mem // min_mem_mb)

    jobs = ssh(remote, "scontrol show job -o", capture=True).stdout.splitlines()
    parsed_jobs = []
    for job in jobs:
        fields = dict(re.findall(r"\b(Comment|JobState|Reason)=(\S+)", job))
        if fields:
            parsed_jobs.append(fields)
    active_ora = sum(
        fields.get("Comment", "").startswith("ora_bp_")
        and fields.get("JobState") in {"RUNNING", "PENDING"}
        for fields in parsed_jobs
    )
    fit_waiting = any(
        fields.get("Comment", "").startswith("fit_bp_")
        and fields.get("JobState") == "PENDING"
        and fields.get("Reason") in {"Resources", "ReqNodeNotAvail"}
        for fields in parsed_jobs
    )
    return free_cpus, free_mem, ora_capacity, active_ora, fit_waiting


def local_ora_slots(args: argparse.Namespace) -> tuple[int, int, int]:
    """Return safe new ORA slots and observed local coverage work."""
    processes = run(["ps", "-eo", "rss=,args="], capture=True).stdout.splitlines()
    fit_processes = [
        command
        for command in processes
        if "--file=scripts/coverage/fit_clampfull.R" in command
    ]
    active_fits = len(fit_processes)
    active_fit_rss_mb = sum(int(command.split(maxsplit=1)[0]) for command in fit_processes) // 1024
    ora_processes = [
        command
        for command in processes
        if "--file=scripts/coverage/run_ora.R" in command
    ]
    active_oras = len(ora_processes)
    active_ora_rss_mb = sum(int(command.split(maxsplit=1)[0]) for command in ora_processes) // 1024
    remaining_mem = max(
        0,
        args.local_total_mem_mb
        - active_fits * args.local_fit_mem_mb
        - active_oras * args.local_ora_mem_mb,
    )
    meminfo = Path("/proc/meminfo").read_text()
    available_match = re.search(r"^MemAvailable:\s+(\d+)\s+kB$", meminfo, re.MULTILINE)
    system_available = (
        int(available_match.group(1)) // 1024
        if available_match
        else args.local_total_mem_mb
    )
    fit_growth_headroom = max(
        0,
        active_fits * args.local_fit_mem_mb - active_fit_rss_mb,
    )
    ora_growth_headroom = max(
        0,
        active_oras * args.local_ora_mem_mb - active_ora_rss_mb,
    )
    usable_system_mem = max(
        0,
        system_available
        - args.local_reserve_mem_mb
        - fit_growth_headroom
        - ora_growth_headroom,
    )
    slots = min(
        args.local_jobs,
        remaining_mem // args.local_ora_mem_mb,
        usable_system_mem // args.local_ora_mem_mb,
    )
    return slots, active_fits, active_oras


def log(message: str) -> None:
    print(f"[route_ora {datetime.now(timezone.utc).isoformat()}] {message}", flush=True)


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temporary.write_text(json.dumps(payload, indent=2) + "\n")
    os.replace(temporary, path)


def discard_unpublished_output(output: Path) -> None:
    """Drop an ORA directory that exists but was never published.

    run_ora.R publishes atomically, so a directory without its `complete` marker
    can only come from an interrupted sync.  Snakemake would treat it as a
    finished `directory()` output and never rebuild it, and run_ora.R refuses to
    overwrite it, so the retry has to start from a clean slate.
    """
    if output.exists() and not (output / "complete").exists():
        log(f"discarding unpublished ORA output: {output}")
        shutil.rmtree(output, ignore_errors=True)


def reclaim_stale_assignments(args: argparse.Namespace) -> None:
    """Make targets abandoned by an earlier router run eligible again.

    A router that is killed mid-target leaves `status: "running"` behind, and the
    pending filter would then skip that target forever.  Nothing else writes
    these manifests, so at startup any `running` entry without a published ORA
    directory belongs to a router that is no longer alive.
    """
    root = Path(args.assignment_root)
    if not root.is_dir():
        return
    for assignment in sorted(root.glob("*.json")):
        try:
            payload = json.loads(assignment.read_text())
        except (OSError, json.JSONDecodeError) as error:
            log(f"unreadable assignment {assignment.name}: {error}")
            continue
        if payload.get("status") != "running":
            continue
        target = Target.from_key(payload.get("target", assignment.stem))
        if target is None:
            log(f"unrecognised assignment target in {assignment.name}")
            continue
        output = Path(args.local_root, target_output(args.ora_root, target))
        if (output / "complete").exists():
            payload.update(status="complete", finished_utc=datetime.now(timezone.utc).isoformat())
            atomic_json(assignment, payload)
            log(f"reclaimed {target.key}: output was already published")
            continue
        discard_unpublished_output(output)
        payload.update(status="interrupted", finished_utc=datetime.now(timezone.utc).isoformat())
        atomic_json(assignment, payload)
        log(f"reclaimed {payload['target']}: retryable again")


def discover_ready(remote: str, remote_root: str, model_root: str) -> list[Target]:
    root = f"{remote_root.rstrip('/')}/{model_root.rstrip('/')}"
    command = f"find {shlex.quote(root)}/archs4 -path '*/seed*/validated.json' -type f -print"
    paths = ssh(remote, command, capture=True).stdout.splitlines()
    targets: list[Target] = []
    for path in paths:
        match = VALIDATED_RE.search(path)
        if match:
            for database in DATABASES:
                targets.append(Target(int(match["fraction"]), int(match["seed"]), database))
    return sorted(targets, key=lambda item: (item.fraction, item.seed, item.database))


def target_output(ora_root: str, target: Target) -> str:
    return f"{ora_root.rstrip('/')}/{target.ora_leaf}"


def sync_model_inputs(args: argparse.Namespace, target: Target) -> None:
    remote_base = (
        f"{args.remote}:{args.remote_root.rstrip('/')}/{args.model_root.rstrip('/')}/"
        f"{target.model_leaf}/"
    )
    local_base = Path(args.local_root, args.model_root, target.model_leaf)
    local_base.mkdir(parents=True, exist_ok=True)
    run([
        "rsync", "-a", "--partial",
        "--include=validated.json", "--include=CLAMPfull_bp/", "--include=CLAMPfull_bp/Z.csv",
        "--include=CLAMPfull_bp/manifest.json", "--exclude=*",
        remote_base, f"{local_base}/",
    ])


def sync_remote_result(args: argparse.Namespace, target: Target) -> None:
    remote_output = (
        f"{args.remote}:{args.remote_root.rstrip('/')}/{target_output(args.ora_root, target)}/"
    )
    local_output = Path(args.local_root, target_output(args.ora_root, target))
    local_output.mkdir(parents=True, exist_ok=True)
    run(["rsync", "-a", "--partial", remote_output, f"{local_output}/"])


def invoke_snakemake(args: argparse.Namespace, target: Target, backend: str) -> None:
    output = target_output(args.ora_root, target)
    if backend == "remote":
        path_prefix = ":".join(shlex.quote(path) for path in args.remote_path_dir)
        environment = f"export PATH={path_prefix}:$PATH; " if path_prefix else ""
        command = (
            environment
            + f"cd {shlex.quote(args.remote_root)} && "
            f"{shlex.quote(args.remote_snakemake)} -s workflow/Snakefile "
            f"--profile {shlex.quote(args.remote_profile)} --nolock {shlex.quote(output)}"
        )
        ssh(args.remote, command)
        sync_remote_result(args, target)
    else:
        sync_model_inputs(args, target)
        run([
            args.local_snakemake,
            "-s", "workflow/Snakefile",
            "--profile", args.local_profile,
            "--nolock", output,
        ], cwd=args.local_root)


def execute_target(args: argparse.Namespace, target: Target, backend: str) -> None:
    assignment = Path(args.assignment_root, f"{target.key}.json")
    now = datetime.now(timezone.utc).isoformat()
    payload = {
        "schema_version": 1,
        "target": target.key,
        "backend": backend,
        "status": "running",
        "assigned_utc": now,
    }
    atomic_json(assignment, payload)
    try:
        invoke_snakemake(args, target, backend)
    except Exception as error:
        discard_unpublished_output(Path(args.local_root, target_output(args.ora_root, target)))
        payload.update(status="failed", error=str(error), finished_utc=datetime.now(timezone.utc).isoformat())
        atomic_json(assignment, payload)
        raise
    payload.update(status="complete", finished_utc=datetime.now(timezone.utc).isoformat())
    atomic_json(assignment, payload)


def choose_backend(args: argparse.Namespace) -> tuple[str, dict]:
    free_cpus, free_mem, ora_capacity, active_ora, fit_waiting = remote_capacity(
        args.remote,
        args.remote_min_free_cpus,
        args.remote_min_free_mem_mb,
    )
    capacity = {
        "free_cpus": free_cpus,
        "free_mem_mb": free_mem,
        "ora_job_capacity": ora_capacity,
        "active_ora_jobs": active_ora,
        "fit_waiting": fit_waiting,
    }
    remote_ok = (
        ora_capacity >= 1
        and active_ora < args.remote_max_ora_jobs
        and not fit_waiting
    )
    return ("remote" if remote_ok else "local"), capacity


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--remote", required=True, help="User-owned SSH host or alias")
    parser.add_argument("--remote-root", required=True)
    parser.add_argument("--local-root", default=".")
    parser.add_argument("--model-root", default="output/03_model_biology/02_archs4/00_coverage/models")
    parser.add_argument("--ora-root", default="output/03_model_biology/02_archs4/00_coverage/ora")
    parser.add_argument("--assignment-root", default="output/03_model_biology/02_archs4/00_coverage/assignments")
    parser.add_argument("--remote-profile", required=True)
    parser.add_argument("--local-profile", default="workflow/profiles/local")
    parser.add_argument("--remote-snakemake", default="snakemake")
    parser.add_argument(
        "--remote-path-dir",
        action="append",
        default=[],
        help="Runtime-only directory prepended to PATH remotely; may be repeated",
    )
    parser.add_argument("--local-snakemake", default="snakemake")
    parser.add_argument("--remote-min-free-cpus", type=int, default=2)
    parser.add_argument("--remote-min-free-mem-mb", type=int, default=16000)
    parser.add_argument("--remote-max-ora-jobs", type=int, default=8)
    parser.add_argument("--local-jobs", type=int, default=8)
    parser.add_argument("--local-total-mem-mb", type=int, default=96000)
    parser.add_argument("--local-fit-mem-mb", type=int, default=32000)
    parser.add_argument("--local-ora-mem-mb", type=int, default=8000)
    parser.add_argument("--local-reserve-mem-mb", type=int, default=16000)
    parser.add_argument("--poll-seconds", type=int, default=300)
    parser.add_argument("--once", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    args.local_root = str(Path(args.local_root).resolve())
    args.assignment_root = str(Path(args.local_root, args.assignment_root))

    if not args.dry_run:
        reclaim_stale_assignments(args)

    while True:
        try:
            ready = discover_ready(args.remote, args.remote_root, args.model_root)
        except Exception as error:
            # Discovery is one ssh call.  A campaign spanning days will lose it
            # occasionally; that is a reason to poll again, not to stop routing.
            if args.dry_run or args.once:
                raise
            log(f"discovery failed, retrying after {args.poll_seconds}s: {error}")
            time.sleep(args.poll_seconds)
            continue

        pending = []
        for target in ready:
            assignment = Path(args.assignment_root, f"{target.key}.json")
            local_complete = Path(args.local_root, target_output(args.ora_root, target), "complete")
            if local_complete.exists():
                continue
            if assignment.exists() and json.loads(assignment.read_text()).get("status") in {"running", "complete"}:
                continue
            pending.append(target)

        if args.dry_run:
            backend, capacity = choose_backend(args)
            local_slots, active_local_fits, active_local_oras = local_ora_slots(args)
            print(json.dumps({
                "capacity": capacity,
                "would_route": backend,
                "local_ora_slots": local_slots,
                "active_local_fits": active_local_fits,
                "active_local_oras": active_local_oras,
                "ready_targets": [item.key for item in pending],
            }, indent=2))
            return

        if pending:
            try:
                preferred, capacity = choose_backend(args)
            except Exception as error:
                if args.once:
                    raise
                log(f"capacity check failed, retrying after {args.poll_seconds}s: {error}")
                time.sleep(args.poll_seconds)
                continue
            remote_slots = 0
            if preferred == "remote":
                remote_slots = min(
                    capacity["ora_job_capacity"],
                    max(0, args.remote_max_ora_jobs - capacity["active_ora_jobs"]),
                )
            local_slots, _, _ = local_ora_slots(args)
            selected = pending[: remote_slots + local_slots]
            workers = min(len(selected), remote_slots + local_slots)
            if workers == 0:
                if args.once:
                    return
                time.sleep(args.poll_seconds)
                continue
            with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
                futures = {}
                for index, target in enumerate(selected):
                    backend = "remote" if index < remote_slots else "local"
                    log(f"routing {target.key} to {backend}")
                    futures[pool.submit(execute_target, args, target, backend)] = target
                for future in concurrent.futures.as_completed(futures):
                    target = futures[future]
                    try:
                        future.result()
                    except Exception as error:
                        # execute_target already recorded the failure in the
                        # assignment manifest, which leaves the target pending
                        # for a later poll.  One bad target must not end the run.
                        log(f"{target.key} failed: {error}")
                    else:
                        log(f"{target.key} complete")

        if args.once:
            return
        time.sleep(args.poll_seconds)


if __name__ == "__main__":
    main()
