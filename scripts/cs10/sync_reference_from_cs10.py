#!/usr/bin/env python3
"""Copy a cs10 reference locally, excluding operational logs, and verify hashes."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import subprocess


EXCLUDED_NAMES = {"run.log", "render.log"}


def _run(args: list[str], *, dry_run: bool) -> None:
    print(" ".join(args))
    if not dry_run:
        subprocess.run(args, check=True)


def _local_hashes(root: Path) -> dict[str, str]:
    return {
        str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(root.rglob("*"))
        if path.is_file() and path.name not in EXCLUDED_NAMES
    }


def _remote_hashes(host: str, remote_dir: str) -> dict[str, str]:
    command = (
        f"cd {remote_dir!s} && find . -type f ! -name run.log ! -name render.log "
        "-print0 | sort -z | xargs -0 sha256sum"
    )
    result = subprocess.run(
        ["ssh", host, command], check=True, text=True, capture_output=True
    )
    return {
        line.split(maxsplit=1)[1].removeprefix("./"): line.split(maxsplit=1)[0]
        for line in result.stdout.splitlines()
    }


def sync(host: str, remote_dir: str, local_dir: Path, *, dry_run: bool = False) -> Path:
    """Synchronize portable reference artifacts via scp without server logs."""
    remote = f"{host}:{remote_dir.rstrip('/')}"
    transfers = [
        ["scp", "-r", f"{remote}/conditions", str(local_dir)],
        ["scp", "-r", f"{remote}/analysis", str(local_dir)],
        *[
            ["scp", f"{remote}/{name}", str(local_dir)]
            for name in (
                "run_manifest.json",
                "reference_manifest.json",
                "summary.csv",
                "campaign_completion.json",
            )
        ],
    ]
    if dry_run:
        for command in transfers:
            _run(command, dry_run=True)
        return local_dir
    local_dir.mkdir(parents=True, exist_ok=False)
    try:
        for command in transfers:
            _run(command, dry_run=False)
        for name in EXCLUDED_NAMES:
            for path in local_dir.rglob(name):
                path.unlink()
        if _local_hashes(local_dir) != _remote_hashes(host, remote_dir):
            raise RuntimeError("cs10 reference checksum mismatch")
    except Exception:
        raise
    return local_dir


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host", required=True)
    parser.add_argument("--remote-dir", required=True)
    parser.add_argument("--local-dir", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)
    print(sync(args.host, args.remote_dir, args.local_dir, dry_run=args.dry_run))


if __name__ == "__main__":
    main()
