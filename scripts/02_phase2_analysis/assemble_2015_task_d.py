#!/usr/bin/env python3
"""Assemble Issue #199 Task D outputs without re-simulating."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.task_d_2015_tau_linked import main


if __name__ == "__main__":
    main()
