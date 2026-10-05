#!/usr/bin/env python3
"""Build a Phase 2 behavior dataset from simulation outputs."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.behavior_dataset import main

if __name__ == "__main__":
    main()
