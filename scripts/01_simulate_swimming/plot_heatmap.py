#!/usr/bin/env python3
"""Plot a configured Phase 2 sweep heatmap."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[2] / "src"))

from sim_swim.analysis.cli_profiles import args_from_profile, load_profile
from sim_swim.analysis.heatmaps import (
    dt_star_torque,
    hook_overstretch,
    local_scale_mode,
    motor_scale_collapse,
)


HEATMAP_MAIN = {
    "motor_scale_collapse": motor_scale_collapse.main,
    "dt_star_torque": dt_star_torque.main,
    "local_scale_mode": local_scale_mode.main,
    "hook_overstretch": hook_overstretch.main,
}


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Heatmap profile YAML under conf/phase2_sweeps/.",
    )
    args, passthrough = parser.parse_known_args(argv)

    profile = load_profile(args.config)
    kind = str(profile.get("kind", "")).strip()
    if kind not in HEATMAP_MAIN:
        choices = ", ".join(sorted(HEATMAP_MAIN))
        raise SystemExit(f"Unknown heatmap kind {kind!r}. Expected one of: {choices}")

    effective_args = args_from_profile(profile) + passthrough
    HEATMAP_MAIN[kind](effective_args)


if __name__ == "__main__":
    main()
