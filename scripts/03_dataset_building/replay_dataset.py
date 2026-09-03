#!/usr/bin/env python3
"""Replay an existing raw run or behavior dataset without re-simulation."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.behavior_dataset_replay import main as replay_dataset
from sim_swim.analysis.phase2_replay import main as replay_run
from sim_swim.analysis.hydrodynamics_replay import main as replay_hydrodynamics


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--dataset-dir", type=Path)
    source.add_argument("--run-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--view", choices=("2d", "3d", "3d+2d"), default="3d+2d")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--flow-overlay", action="store_true")
    args, passthrough = parser.parse_known_args(argv)
    if args.flow_overlay:
        if args.run_dir is None:
            parser.error("--flow-overlay requires --run-dir")
        replay_args = ["--run-dir", str(args.run_dir), *passthrough]
        if args.output_dir is not None:
            replay_args.extend(["--output-dir", str(args.output_dir)])
        replay_hydrodynamics(replay_args)
        return
    if args.dataset_dir is not None:
        replay_args = ["--dataset-dir", str(args.dataset_dir)]
        if args.view == "2d":
            replay_args.append("--no-render-3d")
        elif args.view == "3d":
            replay_args.append("--no-render-2d")
        if args.output_dir is not None:
            replay_args.extend(["--output-dir", str(args.output_dir)])
        replay_args.extend(passthrough)
        replay_dataset(replay_args)
        return
    replay_args = [
        "--run-dir",
        str(args.run_dir),
        "--mode",
        "render-only",
        "--view",
        args.view,
    ]
    if args.output_dir is not None:
        replay_args.extend(["--output-dir", str(args.output_dir)])
    if args.overwrite:
        replay_args.append("--overwrite")
    replay_args.extend(passthrough)
    replay_run(replay_args)


if __name__ == "__main__":
    main()
