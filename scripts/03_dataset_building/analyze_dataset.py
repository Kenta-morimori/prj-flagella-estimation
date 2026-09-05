#!/usr/bin/env python3
"""Create reusable 2D and 3D analyses for an existing behavior dataset."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.behavior_dataset_distributions import (
    analyze_dataset as analyze_3d,
)
from sim_swim.analysis.behavior_dataset_separability import analyze_2d_separability
from sim_swim.analysis import common_analysis

HEATMAP_MAIN = common_analysis.HEATMAP_MAIN
dispatch = common_analysis.dispatch


def main(argv: list[str] | None = None) -> None:
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    if "--analysis-kind" in raw_argv:
        index = raw_argv.index("--analysis-kind")
        try:
            analysis_kind = raw_argv[index + 1]
        except IndexError as error:
            raise SystemExit("--analysis-kind requires a value") from error
        dispatch(analysis_kind, raw_argv[:index] + raw_argv[index + 2 :])
        return
    if raw_argv and (
        any(value.startswith("config=") for value in raw_argv)
        or "--config" in raw_argv
        or any(
            value in {"list_profiles=true", "list_canonical_profiles=true"}
            for value in raw_argv
        )
    ):
        dispatch("heatmap", raw_argv)
        return
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--view", choices=("2d", "3d", "both"), default="both")
    parser.add_argument("--n-flagella", default="1,2,3")
    parser.add_argument("--include-non-ml-candidates", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument(
        "--analysis-kind",
        choices=(
            "heatmap",
            "2010-torque-dt",
            "issue61-2015-1tau",
            "2010-fixed-performance",
            "2015-stage-a",
            "spring-formulations",
            "task-d-2015",
            "diagnose-158",
            "partial-generic",
            "plan-2010-torque-dt",
            "plan-reference-torque",
            "visualize-2010-torque-dt",
            "phase-seed-2d",
            "distributions",
            "run-summary",
        ),
        help="Run a reusable raw-run or campaign analysis instead of dataset analysis.",
    )
    args = parser.parse_args(raw_argv)
    dataset_dir = args.dataset_dir.resolve()
    output_dir = args.output_dir or dataset_dir / "analysis" / "common"
    outputs: dict[str, object] = {}
    if args.view in {"3d", "both"}:
        # The distribution workflow reads the source 3D summary features recorded
        # in the behavior dataset and writes its own reproducible subdirectory.
        outputs["3d"] = {
            key: [str(path) for path in value]
            if isinstance(value, list)
            else str(value)
            for key, value in analyze_3d(
                dataset_dir=dataset_dir, overwrite=args.overwrite
            ).items()
        }
    output_dir.mkdir(parents=True, exist_ok=True)
    if args.view in {"2d", "both"}:
        n_flagella = tuple(int(value) for value in args.n_flagella.split(",") if value)
        outputs["2d"] = str(
            analyze_2d_separability(
                dataset_dir=dataset_dir,
                output_dir=output_dir / "projection_2d",
                n_flagella=n_flagella,
                ml_candidates_only=not args.include_non_ml_candidates,
                overwrite=args.overwrite,
            )
        )
    (output_dir / "manifest.json").write_text(
        json.dumps(
            {
                "kind": "common_dataset_analysis",
                "dataset_dir": str(dataset_dir),
                "view": args.view,
                "outputs": outputs,
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(output_dir)


if __name__ == "__main__":
    main()
