#!/usr/bin/env python3
"""Compare compact Phase 2 serial or parallel qualification outputs."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[2] / "src"))

from sim_swim.analysis.qualification import (
    DEFAULT_ATOL,
    DEFAULT_RTOL,
    compare_campaigns,
    compare_parallel_job,
    write_report,
)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--left", type=Path, help="First campaign directory or run_manifest.json"
    )
    parser.add_argument(
        "--right", type=Path, help="Second campaign directory or run_manifest.json"
    )
    parser.add_argument(
        "--parallel-job-manifest",
        type=Path,
        help="job_manifest.json for serial-versus-parallel comparison",
    )
    parser.add_argument(
        "--serial",
        action="append",
        default=[],
        metavar="CONFIG=PATH",
        help="Serial reference, repeat once for each parallel config",
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--atol", type=float, default=DEFAULT_ATOL)
    parser.add_argument("--rtol", type=float, default=DEFAULT_RTOL)
    args = parser.parse_args(argv)
    if args.parallel_job_manifest:
        if args.left or args.right:
            parser.error(
                "--parallel-job-manifest cannot be combined with --left/--right"
            )
        serial_campaigns: dict[str, Path] = {}
        for item in args.serial:
            config, separator, raw_path = item.partition("=")
            if not separator or not config or not raw_path:
                parser.error("--serial must use CONFIG=PATH")
            serial_campaigns[config] = Path(raw_path)
        report = compare_parallel_job(
            args.parallel_job_manifest, serial_campaigns, atol=args.atol, rtol=args.rtol
        )
    else:
        if not args.left or not args.right:
            parser.error(
                "--left and --right are required without --parallel-job-manifest"
            )
        report = compare_campaigns(
            args.left, args.right, atol=args.atol, rtol=args.rtol
        )
    json_path, csv_path = write_report(report, args.output_dir)
    print(report["status"])
    print(json_path)
    print(csv_path)
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
