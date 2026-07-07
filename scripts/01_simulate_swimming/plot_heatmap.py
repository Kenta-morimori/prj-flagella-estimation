#!/usr/bin/env python3
"""Plot a configured Phase 2 sweep heatmap."""

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
)
from sim_swim.analysis.heatmaps import (
    dt_star_torque,
    hook_overstretch,
    local_scale_mode,
    motor_scale_collapse,
    shape_stability_grid,
)


HEATMAP_MAIN = {
    "motor_scale_collapse": motor_scale_collapse.main,
    "dt_star_torque": dt_star_torque.main,
    "local_scale_mode": local_scale_mode.main,
    "shape_stability_grid": shape_stability_grid.main,
    "hook_overstretch": hook_overstretch.main,
}


def _has_option(args: list[str], option_name: str) -> bool:
    return any(arg == option_name or arg.startswith(f"{option_name}=") for arg in args)


def _first_option_value(args: list[str], option_name: str) -> str | None:
    for index, arg in enumerate(args):
        if arg == option_name and index + 1 < len(args):
            return args[index + 1]
        if arg.startswith(f"{option_name}="):
            return arg.split("=", 1)[1]
    return None


def _with_default_output_dir(args: list[str]) -> list[str]:
    if _has_option(args, "--output-dir"):
        return args
    summary_csv = _first_option_value(args, "--summary-csv")
    if summary_csv is None:
        return args
    return [*args, "--output-dir", str(Path(summary_csv).parent / "plots")]


def main(argv: list[str] | None = None) -> None:
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    config_from_key, parser_argv = split_config_key(raw_argv)

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Heatmap profile YAML under conf/phase2_sweeps/.",
    )
    args, passthrough = parser.parse_known_args(parser_argv)
    if config_from_key is not None and args.config is not None:
        parser.error("Use either config=PATH or --config PATH (not both)")
    config = config_from_key or args.config
    if config is None:
        parser.error("config=PATH or --config PATH is required")

    profile = load_profile(config)
    kind = str(profile.get("kind", "")).strip()
    if kind not in HEATMAP_MAIN:
        choices = ", ".join(sorted(HEATMAP_MAIN))
        raise SystemExit(f"Unknown heatmap kind {kind!r}. Expected one of: {choices}")

    effective_args = _with_default_output_dir(
        args_from_profile(profile) + key_value_args_to_cli_args(passthrough)
    )
    HEATMAP_MAIN[kind](effective_args)


if __name__ == "__main__":
    main()
