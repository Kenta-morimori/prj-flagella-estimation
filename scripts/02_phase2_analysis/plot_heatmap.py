#!/usr/bin/env python3
"""Plot a configured Phase 2 sweep heatmap."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.phase2_heatmap import HEATMAP_MAIN, main

__all__ = ["HEATMAP_MAIN", "main"]


if __name__ == "__main__":
    main()
