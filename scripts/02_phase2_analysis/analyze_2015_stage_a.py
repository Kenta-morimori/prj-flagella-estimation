#!/usr/bin/env python3
"""Analyze Issue #168 motor-off pilot and optional motor-on Stage A results."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.stage_a_2015_workflow import main


if __name__ == "__main__":
    main()
