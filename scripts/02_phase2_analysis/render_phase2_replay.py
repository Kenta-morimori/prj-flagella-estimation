#!/usr/bin/env python3
"""Render or plot Phase 2 replay diagnostics from existing simulation outputs."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.phase2_replay import main


if __name__ == "__main__":
    main()
