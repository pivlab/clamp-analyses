#!/usr/bin/env python3

# Publishes a full-data GO:BP model through links and a manifest

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

PAYLOAD = ("CLAMPfull_bp.rds", "B.csv", "Z.csv", "summary.csv")


# Calculates a file checksum.
def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


# Creates a hard link or symbolic link for a model artifact.
def link(src: Path, dst: Path) -> str:
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    try:
        os.link(src, dst)
        return "hardlink"
    except OSError:
        dst.symlink_to(src.resolve())
        return "symlink"


# Publishes the requested model directory.
def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--model-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--validated", type=Path, default=None)
    args = parser.parse_args()

    model_dir: Path = args.model_dir
    missing = [name for name in PAYLOAD if not (model_dir / name).is_file()]
    if missing:
        raise SystemExit(f"{args.dataset}: missing model artifacts: {', '.join(missing)}")

    source_manifest = model_dir / "manifest.json"
    if not source_manifest.is_file():
        raise SystemExit(f"{args.dataset}: no manifest.json in {model_dir}")
    manifest = json.loads(source_manifest.read_text())
    if manifest.get("status") != "complete":
        raise SystemExit(f"{args.dataset}: source manifest is not complete")

    payload_dir = args.out_dir / "CLAMPfull_bp"
    payload_dir.mkdir(parents=True, exist_ok=True)

    modes = {}
    for name in PAYLOAD:
        modes[name] = link(model_dir / name, payload_dir / name)

    modes["CLAMPfull_bp.rds@top"] = link(
        model_dir / "CLAMPfull_bp.rds", args.out_dir / "CLAMPfull_bp.rds"
    )

    published = {
        "schema_version": 1,
        "dataset": args.dataset,
        "model": "CLAMPfull_bp",
        "prior": manifest.get("prior"),
        "genes": manifest.get("genes"),
        "samples": manifest.get("samples"),
        "latent_variables": manifest.get("latent_variables"),
        "rng_seed": manifest.get("rng_seed"),
        "source": {
            "model_dir": str(model_dir.resolve()),
            "fraction": manifest.get("fraction"),
            "seed_index": manifest.get("seed_index"),
            "manifest_sha256": sha256(source_manifest),
        },
        "z_sha256": sha256(model_dir / "Z.csv"),
        "link_mode": modes,
    }
    if args.validated is not None and args.validated.is_file():
        published["validated"] = json.loads(args.validated.read_text()).get("status")

    (args.out_dir / "manifest.json").write_text(json.dumps(published, indent=2) + "\n")
    print(f"Published {args.dataset} -> {args.out_dir} ({set(modes.values())})")


if __name__ == "__main__":
    main()
