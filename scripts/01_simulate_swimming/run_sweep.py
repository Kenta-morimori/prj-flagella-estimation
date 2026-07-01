#!/usr/bin/env python3
"""Run a configured Phase 2 simulation sweep."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[2] / "src"))

from sim_swim.analysis.cli_profiles import (
    args_from_profile,
    key_value_args_to_cli_args,
    load_profile,
    split_config_key,
    sweep_aliases,
)
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
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    config_from_key, parser_argv = split_config_key(raw_argv)
    parser_argv = key_value_args_to_cli_args(parser_argv)

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Sweep profile YAML under conf/phase2_sweeps/.",
    )
    parser.add_argument(
        "--list-kind",
        action="store_true",
        help="Print the profile kind and exit without running conditions.",
    )
    args, passthrough = parser.parse_known_args(parser_argv)
    if config_from_key is not None and args.config is not None:
        parser.error("Use either config=PATH or --config PATH (not both)")
    config = config_from_key or args.config
    if config is None:
        parser.error("config=PATH or --config PATH is required")

    profile = load_profile(config)
    kind = str(profile.get("kind", "")).strip()
    if kind not in SWEEP_MAIN:
        choices = ", ".join(sorted(SWEEP_MAIN))
        raise SystemExit(f"Unknown sweep kind {kind!r}. Expected one of: {choices}")
    if args.list_kind:
        print(kind)
        return

    effective_args = args_from_profile(profile) + key_value_args_to_cli_args(
        passthrough,
        aliases=sweep_aliases(kind),
    )
    SWEEP_MAIN[kind](effective_args)


if __name__ == "__main__":
    main()
