#!/usr/bin/env python3
"""Run a configured Phase 2 simulation sweep."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[2] / "src"))

from sim_swim.analysis.cli_profiles import (
    args_from_profile,
    format_profile_description,
    format_profile_listing,
    key_value_args_to_cli_args,
    list_profile_entries,
    load_profile_entry,
    split_config_key,
    sweep_aliases,
    validate_profile_role,
)
from sim_swim.analysis.sweeps import (
    bundling_alignment,
    generic_multi_run,
    hook_overstretch,
    motor_scale,
    shape_stability_grid,
    single_flagellum_torque,
)


SWEEP_MAIN = {
    "motor_scale": motor_scale.main,
    "single_flagellum_torque": single_flagellum_torque.main,
    "bundling_alignment": bundling_alignment.main,
    "shape_stability_grid": shape_stability_grid.main,
    "hook_overstretch": hook_overstretch.main,
    "generic_multi_run": generic_multi_run.main,
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
    parser.add_argument(
        "--list-profiles",
        action="store_true",
        help="List available sweep profiles and exit.",
    )
    parser.add_argument(
        "--list-canonical-profiles",
        action="store_true",
        help="List canonical sweep profiles and exit.",
    )
    parser.add_argument(
        "--describe-profile",
        action="store_true",
        help="Print profile metadata for the selected config and exit.",
    )
    args, passthrough = parser.parse_known_args(parser_argv)
    if config_from_key is not None and args.config is not None:
        parser.error("Use either config=PATH or --config PATH (not both)")
    if args.list_profiles or args.list_canonical_profiles:
        entries = list_profile_entries(
            role="sweep", canonical_only=args.list_canonical_profiles
        )
        for line in format_profile_listing(entries):
            print(line)
        return
    config = config_from_key or args.config
    if config is None:
        parser.error("config=PATH or --config PATH is required")

    entry = load_profile_entry(config)
    kind = entry["kind"]
    if kind not in SWEEP_MAIN:
        choices = ", ".join(sorted(SWEEP_MAIN))
        raise SystemExit(f"Unknown sweep kind {kind!r}. Expected one of: {choices}")
    if args.describe_profile:
        for line in format_profile_description(entry, list_profile_entries()):
            print(line)
        return
    if args.list_kind:
        print(kind)
        return
    validate_profile_role(entry, "sweep")

    if kind == "generic_multi_run":
        effective_args = [
            "--campaign-config",
            str(config),
        ] + key_value_args_to_cli_args(passthrough)
    else:
        effective_args = args_from_profile(entry) + key_value_args_to_cli_args(
            passthrough,
            aliases=sweep_aliases(kind),
        )
    SWEEP_MAIN[kind](effective_args)


if __name__ == "__main__":
    main()
