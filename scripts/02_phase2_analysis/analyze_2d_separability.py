#!/usr/bin/env python3
"""Analyze 2D separability for an existing Phase 2 behavior dataset."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.behavior_dataset_separability import main


if __name__ == "__main__":
    main()
