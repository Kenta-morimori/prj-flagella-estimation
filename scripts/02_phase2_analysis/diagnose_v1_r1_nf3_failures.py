#!/usr/bin/env python3
"""Diagnose Phase 2 Issue #158 v1 r1 n=3 long-run failures."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.phase2_158_diagnostics import (
    DEFAULT_FAIL_CONDITION_IDS,
    Phase2158DiagnosticConfig,
    analyze_phase2_158_diagnostics,
    default_output_dir,
)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--campaign-config",
        type=Path,
        default=Path("conf/phase2_multi_run/flagella_count_duration_3s_r1.yaml"),
        help="Phase 2 multi-run campaign config used to reconstruct condition topology.",
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=Path("outputs/phase2_multi_run/flagella_count_duration_3s_r1"),
        help="Existing raw multi-run output directory.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory. Defaults to outputs/YYYY-MM-DD/HHMMSS/...",
    )
    parser.add_argument(
        "--fail-condition-id",
        action="append",
        default=None,
        help="Condition ID to include in first-fail event plots. Can be repeated.",
    )
    parser.add_argument(
        "--first-fail-window-s",
        type=float,
        default=0.25,
        help="Seconds before/after first fail to include in event tables and plots.",
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "overrides",
        nargs="*",
        help="Optional KEY=VALUE campaign overrides for analysis reproducibility.",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    output_dir = (
        args.output_dir if args.output_dir is not None else default_output_dir()
    )
    result = analyze_phase2_158_diagnostics(
        Phase2158DiagnosticConfig(
            campaign_config=args.campaign_config,
            input_dir=args.input_dir,
            output_dir=output_dir,
            fail_condition_ids=tuple(
                args.fail_condition_id or DEFAULT_FAIL_CONDITION_IDS
            ),
            first_fail_window_s=float(args.first_fail_window_s),
            overwrite=bool(args.overwrite),
            cli_overrides=tuple(args.overrides),
        )
    )
    print(result)


if __name__ == "__main__":
    main()
