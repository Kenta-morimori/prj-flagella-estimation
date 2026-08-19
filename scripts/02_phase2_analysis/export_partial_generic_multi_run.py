#!/usr/bin/env python3
"""Export completed conditions from an interrupted generic multi-run campaign."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.partial_generic_multi_run import main


if __name__ == "__main__":
    main()
