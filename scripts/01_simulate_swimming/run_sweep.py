#!/usr/bin/env python3
"""Run a configured Phase 2 simulation sweep."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[2] / "src"))

from sim_swim.analysis.cli_profiles import args_from_profile, load_profile
from sim_swim.analysis.sweeps import (
    bundling_alignment,
    hook_overstretch,
    motor_scale,
    single_flagellum_torque,
)


SWEEP_MAIN = {
    "motor_scale": motor_scale.main,
    "single_flagellum_torque": single_flagellum_torque.main,
    "bundling_alignment": bundling_alignment.main,
    "hook_overstretch": hook_overstretch.main,
}


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        required=True,
        help="Sweep profile YAML under conf/phase2_sweeps/.",
    )
    parser.add_argument(
        "--list-kind",
        action="store_true",
        help="Print the profile kind and exit without running conditions.",
    )
    args, passthrough = parser.parse_known_args(argv)

    profile = load_profile(args.config)
    kind = str(profile.get("kind", "")).strip()
    if kind not in SWEEP_MAIN:
        choices = ", ".join(sorted(SWEEP_MAIN))
        raise SystemExit(f"Unknown sweep kind {kind!r}. Expected one of: {choices}")
    if args.list_kind:
        print(kind)
        return

    effective_args = args_from_profile(profile) + passthrough
    SWEEP_MAIN[kind](effective_args)


if __name__ == "__main__":
    main()
