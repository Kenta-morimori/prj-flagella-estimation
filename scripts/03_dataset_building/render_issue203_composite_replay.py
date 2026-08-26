#!/usr/bin/env python3
"""Render one Mac-side #203 composite replay without re-simulation."""

from __future__ import annotations
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from sim_swim.analysis.issue203_composite_replay import main

if __name__ == "__main__":
    main()
