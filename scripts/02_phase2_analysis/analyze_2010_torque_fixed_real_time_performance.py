#!/usr/bin/env python3
"""Rebuild Issue #193 fixed-real-time performance outputs without simulation."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[2] / "src"))

from sim_swim.analysis.torque_dt_stability_campaign import (
    render_fixed_real_time_qualitative_replay,
    summarize_campaign,
)


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
    parser.add_argument(
        "--allow-incomplete",
        action="store_true",
        help=(
            "Record uncompleted conditions as exclusions. Required for a "
            "partial qualitative replay; never enables a performance conclusion."
        ),
    )
    parser.add_argument(
        "--render-qualitative-replay",
        action="store_true",
        help="Render only completed conditions from existing archives.",
    )
    parser.add_argument("--fps-out-3d", type=float, default=20.0)
    parser.add_argument("--target-frame-count", type=int, default=101)
    args = parser.parse_args(argv)
    summarize_campaign(
        args.run_dir,
        config_path=args.config,
        allow_incomplete=args.allow_incomplete,
    )
    if args.render_qualitative_replay:
        output_dir = render_fixed_real_time_qualitative_replay(
            args.run_dir,
            config_path=args.config,
            allow_incomplete=args.allow_incomplete,
            fps_out_3d=args.fps_out_3d,
            target_frame_count=args.target_frame_count,
        )
        print(output_dir)
    print(args.run_dir / "performance_summary.csv")


if __name__ == "__main__":
    main()
