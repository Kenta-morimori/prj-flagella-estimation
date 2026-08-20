#!/usr/bin/env python3
"""Reject active references to scripts deleted by the current change set."""

from __future__ import annotations

import argparse
import subprocess


HISTORICAL_PREFIX = "docs/codex-runs/"


def _git(*args: str) -> str:
    return subprocess.check_output(["git", *args], text=True).strip()


def deleted_scripts(base: str) -> list[str]:
    rows = _git("diff", "--name-status", f"{base}...HEAD").splitlines()
    return [
        row.split("\t", 1)[1]
        for row in rows
        if row.startswith("D\t") and row.split("\t", 1)[1].startswith("scripts/")
    ]


def active_references(path: str) -> list[str]:
    result = subprocess.run(
        ["git", "grep", "-n", "--", path], text=True, capture_output=True, check=False
    )
    if result.returncode not in {0, 1}:
        raise subprocess.CalledProcessError(
            result.returncode, result.args, output=result.stdout, stderr=result.stderr
        )
    matches = result.stdout.strip().splitlines() if result.stdout else []
    return [row for row in matches if not row.startswith(HISTORICAL_PREFIX)]


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", default="origin/main")
    args = parser.parse_args(argv)
    violations = {path: active_references(path) for path in deleted_scripts(args.base)}
    violations = {path: rows for path, rows in violations.items() if rows}
    if violations:
        details = "\n".join(
            f"{path}:\n" + "\n".join(f"  {row}" for row in rows)
            for path, rows in violations.items()
        )
        raise SystemExit(
            "Deleted issue-specific scripts still have active consumers. "
            "Migrate them before deletion:\n" + details
        )


if __name__ == "__main__":
    main()
