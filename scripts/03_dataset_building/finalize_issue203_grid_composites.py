#!/usr/bin/env python3
"""Retain only the verified n-flagella grid composites for Issue #203."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.issue203_artifact_layout import finalize_grid_composites


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--uniform-root", type=Path, required=True)
    parser.add_argument("--diffusive-root", type=Path, required=True)
    parser.add_argument("--apply", action="store_true")
    args = parser.parse_args()
    print(
        json.dumps(
            {
                "uniform_grid_count": finalize_grid_composites(
                    args.uniform_root, "uniform", apply=args.apply
                ),
                "diffusive_grid_count": finalize_grid_composites(
                    args.diffusive_root, "diffusive", apply=args.apply
                ),
                "applied": args.apply,
            },
            ensure_ascii=False,
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
