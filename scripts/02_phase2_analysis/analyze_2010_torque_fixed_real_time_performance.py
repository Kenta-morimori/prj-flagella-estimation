#!/usr/bin/env python3
"""Rebuild Issue #193 fixed-real-time performance outputs without simulation."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[2] / "src"))

from sim_swim.analysis.torque_dt_stability_campaign import summarize_campaign


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument(
        "--config",
        type=Path,
        default=Path(
            "conf/phase2_multi_run/"
            "2010_project_torque_fixed_real_time_0p5s_performance.yaml"
        ),
    )
    args = parser.parse_args(argv)
    summarize_campaign(args.run_dir, config_path=args.config)
    print(args.run_dir / "performance_summary.csv")


if __name__ == "__main__":
    main()
