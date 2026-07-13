#!/usr/bin/env python3
"""Run a configured Phase 2 generic multi-run campaign."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[2] / "src"))

from sim_swim.analysis.cli_profiles import (
    format_profile_description,
    format_profile_listing,
    key_value_args_to_cli_args,
    list_profile_entries,
    load_profile_entry,
    split_config_key,
)
from sim_swim.analysis.sweeps import generic_multi_run


def _campaign_entries(*, canonical_only: bool = False) -> list[dict[str, object]]:
    return [
        entry
        for entry in list_profile_entries(role="sweep", canonical_only=canonical_only)
        if entry["kind"] == "generic_multi_run"
    ]


def main(argv: list[str] | None = None) -> None:
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    config_from_key, parser_argv = split_config_key(raw_argv)
    parser_argv = key_value_args_to_cli_args(parser_argv)

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Generic multi-run profile YAML under conf/phase2_multi_run/.",
    )
    parser.add_argument(
        "--list-kind",
        action="store_true",
        help="Print the profile kind and exit without running conditions.",
    )
    parser.add_argument(
        "--list-profiles",
        action="store_true",
        help="List available generic multi-run profiles and exit.",
    )
    parser.add_argument(
        "--list-canonical-profiles",
        action="store_true",
        help="List canonical generic multi-run profiles and exit.",
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
        entries = _campaign_entries(canonical_only=args.list_canonical_profiles)
        for line in format_profile_listing(entries):
            print(line)
        return

    config = config_from_key or args.config
    if config is None:
        parser.error("config=PATH or --config PATH is required")

    entry = load_profile_entry(config)
    if entry["kind"] != "generic_multi_run":
        raise SystemExit(
            f"Profile {entry['path']!r} is kind {entry['kind']!r}; "
            "use run_sweep.py for task-specific sweeps."
        )
    if args.describe_profile:
        for line in format_profile_description(entry, _campaign_entries()):
            print(line)
        return
    if args.list_kind:
        print(entry["kind"])
        return

    effective_args = ["--campaign-config", str(config)] + key_value_args_to_cli_args(
        passthrough
    )
    generic_multi_run.main(effective_args)


if __name__ == "__main__":
    main()
