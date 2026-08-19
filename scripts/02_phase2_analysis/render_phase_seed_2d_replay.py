#!/usr/bin/env python3
"""Render a full-duration body-only 2D grid for phase-seed comparison."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.phase_seed_2d_replay import main


if __name__ == "__main__":
    main()
