#!/usr/bin/env python3
"""Analyze reusable 3D/2D motion features from a generic multi-run campaign."""

from __future__ import annotations
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from sim_swim.analysis.motion_feature_study import (
    analyze_motion_feature_study,
    load_config,
)


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("overrides", nargs="*")
    args = parser.parse_args(argv)
    print(analyze_motion_feature_study(load_config(args.config, args.overrides)))


if __name__ == "__main__":
    main()
