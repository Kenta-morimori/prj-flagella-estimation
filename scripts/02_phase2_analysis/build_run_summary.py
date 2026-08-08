#!/usr/bin/env python3
"""既存 Phase 2 diagnostics から軽量な run_summary.json を生成する。"""

from __future__ import annotations

import argparse
from pathlib import Path

from sim_swim.analysis.run_summary import write_run_summary


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    print(write_run_summary(args.input_dir, overwrite=args.overwrite))


if __name__ == "__main__":
    main()
