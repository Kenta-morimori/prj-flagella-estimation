#!/usr/bin/env python3
"""Run independent Phase 2 sweep profiles in bounded parallel subprocesses."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[2] / "src"))

from sim_swim.analysis.cli_profiles import key_value_args_to_cli_args, split_config_key
from sim_swim.analysis.parallel_job import (
    build_plan,
    job_output_root,
    load_parallel_job,
    resolve_execution,
    run_parallel_job,
)


def _parse_workers(value: str) -> int | str:
    if value == "auto":
        return value
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            "max_workers must be a positive integer or 'auto'"
        ) from exc
    if parsed < 1:
        raise argparse.ArgumentTypeError("max_workers must be positive")
    return parsed


def main(argv: list[str] | None = None) -> int:
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    config_from_key, parser_argv = split_config_key(raw_argv)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--max-workers", type=_parse_workers, default=None)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(key_value_args_to_cli_args(parser_argv))
    if config_from_key is not None and args.config is not None:
        parser.error("Use either config=PATH or --config PATH (not both)")
    if (config := config_from_key or args.config) is None:
        parser.error("config=PATH or --config PATH is required")
    job = load_parallel_job(config)
    execution = resolve_execution(job, args.max_workers)
    if args.dry_run:
        print(
            json.dumps(
                build_plan(job, execution, job_output_root(job)),
                ensure_ascii=False,
                indent=2,
            )
        )
        return 0
    manifest = run_parallel_job(job, execution)
    print(manifest["output_root"])
    return 0 if manifest["status"] == "succeeded" else 1


if __name__ == "__main__":
    raise SystemExit(main())
