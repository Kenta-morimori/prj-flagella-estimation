"""Estimate a longer run from a completed parallel generic campaign."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "src"))

from sim_swim.analysis.runtime_projection import (  # noqa: E402
    estimate_parallel_runtime,
    write_runtime_projection,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--job-root", type=Path, required=True)
    parser.add_argument("--target-duration-s", type=float, required=True)
    parser.add_argument(
        "--conditions", required=True, help="Comma-separated condition IDs"
    )
    parser.add_argument(
        "--historical-performance", action="append", default=[], metavar="ID=PATH"
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    historical = {}
    for item in args.historical_performance:
        condition_id, separator, path = item.partition("=")
        if not separator or not condition_id or not path or condition_id in historical:
            parser.error(f"invalid --historical-performance: {item}")
        historical[condition_id] = Path(path)
    conditions = tuple(args.conditions.split(","))
    if (
        not conditions
        or any(not item for item in conditions)
        or len(set(conditions)) != len(conditions)
    ):
        parser.error("--conditions must list distinct nonempty IDs")
    if set(historical) - set(conditions):
        parser.error("historical ID is not in --conditions")
    result = estimate_parallel_runtime(
        args.job_root,
        target_duration_s=args.target_duration_s,
        expected_conditions=conditions,
        historical_performance=historical,
    )
    write_runtime_projection(result, args.output_dir)
    print(args.output_dir / "runtime_projection.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
