#!/usr/bin/env python3
"""Rebuild Issue #61 2010 torque-linked dt QC from an existing campaign."""

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
            "conf/phase2_multi_run/2010_project_torque_dt_initial_screen.yaml"
        ),
    )
    args = parser.parse_args(argv)
    summarize_campaign(args.run_dir, config_path=args.config)
    print(args.run_dir / "qc_summary.json")


if __name__ == "__main__":
    main()
