#!/usr/bin/env python3
"""Compare the 2 s and 5 s Issue #215 motion-feature references."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from sim_swim.analysis.axis_angle_comparison import main

if __name__ == "__main__":
    main()
