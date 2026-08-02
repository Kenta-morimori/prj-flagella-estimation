#!/usr/bin/env python3
"""Generate Issue #163 spring formulation comparison and default decision."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[2] / "src"))

from sim_swim.analysis.cli_profiles import key_value_args_to_cli_args
from sim_swim.analysis.spring_formulation_comparison import (
    write_decision_artifacts,
    write_force_extension_artifacts,
)
from sim_swim.core.run_context import init_run


def main(argv: list[str] | None = None) -> Path:
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--motor-off-summary", type=Path, required=True)
    parser.add_argument("--motor-on-summary", type=Path, required=True)
    parser.add_argument(
        "--output-base-dir",
        type=Path,
        default=Path("outputs/phase2_potential_comparison"),
    )
    parser.add_argument("--samples", type=int, default=199)
    args = parser.parse_args(key_value_args_to_cli_args(raw_argv))

    ctx = init_run(
        args.output_base_dir,
        input_info={
            "motor_off_summary": str(args.motor_off_summary),
            "motor_on_summary": str(args.motor_on_summary),
            "samples": args.samples,
        },
    )
    csv_path, png_path = write_force_extension_artifacts(
        ctx.out.root,
        samples=args.samples,
    )
    json_path, markdown_path, decision = write_decision_artifacts(
        ctx.out.root,
        motor_off_summary=args.motor_off_summary,
        motor_on_summary=args.motor_on_summary,
    )
    manifest_path = ctx.out.root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["kind"] = "spring_formulation_comparison"
    manifest["outputs"].update(
        {
            "force_extension_csv": str(csv_path),
            "force_extension_png": str(png_path),
            "default_decision_json": str(json_path),
            "default_decision_markdown": str(markdown_path),
        }
    )
    manifest["decision"] = {
        "selected_default": decision["selected_default"],
        "reason": decision["reason"],
    }
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    ctx.logger.info("Selected default: %s", decision["selected_default"])
    ctx.logger.info("Wrote comparison artifacts: %s", ctx.out.root)
    print(ctx.out.root)
    return ctx.out.root


if __name__ == "__main__":
    main()
