#!/usr/bin/env python3
"""Analyze Issue #168 motor-off pilot and optional motor-on Stage A results."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[2] / "src"))

from sim_swim.analysis.cli_profiles import key_value_args_to_cli_args
from sim_swim.analysis.stage_a_2015_analysis import write_analysis
from sim_swim.core.run_context import init_run


def main(argv: list[str] | None = None) -> Path:
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--motor-off-run", type=Path, required=True)
    parser.add_argument("--motor-on-run", type=Path, default=None)
    parser.add_argument("--threshold-contract", type=Path, default=None)
    parser.add_argument("--coarse-motor-off-run", type=Path, default=None)
    parser.add_argument("--coarse-motor-on-run", type=Path, default=None)
    parser.add_argument("--fine-motor-on-short-run", type=Path, default=None)
    parser.add_argument("--output-base-dir", type=Path, default=Path("outputs"))
    args = parser.parse_args(key_value_args_to_cli_args(raw_argv))

    ctx = init_run(
        args.output_base_dir,
        input_info={
            "kind": "stage_a_2015_analysis",
            "motor_off_run": str(args.motor_off_run),
            "motor_on_run": str(args.motor_on_run) if args.motor_on_run else None,
            "threshold_contract": (
                str(args.threshold_contract) if args.threshold_contract else None
            ),
            "coarse_motor_off_run": (
                str(args.coarse_motor_off_run) if args.coarse_motor_off_run else None
            ),
            "coarse_motor_on_run": (
                str(args.coarse_motor_on_run) if args.coarse_motor_on_run else None
            ),
            "fine_motor_on_short_run": (
                str(args.fine_motor_on_short_run)
                if args.fine_motor_on_short_run
                else None
            ),
        },
    )
    outputs = write_analysis(
        ctx.out.root,
        motor_off_run=args.motor_off_run,
        motor_on_run=args.motor_on_run,
        threshold_contract=args.threshold_contract,
        coarse_motor_off_run=args.coarse_motor_off_run,
        coarse_motor_on_run=args.coarse_motor_on_run,
        fine_motor_on_short_run=args.fine_motor_on_short_run,
    )
    manifest_path = ctx.out.root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["kind"] = "stage_a_2015_analysis"
    manifest["outputs"].update({key: str(value) for key, value in outputs.items()})
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    ctx.logger.info("Wrote Stage A analysis: %s", ctx.out.root)
    print(ctx.out.root)
    return ctx.out.root


if __name__ == "__main__":
    main()
