#!/usr/bin/env python3
"""Create Issue #61 heatmaps from an existing torque-dt campaign."""

from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[2] / "src"))

from sim_swim.analysis.torque_dt_stability_visuals import main


if __name__ == "__main__":
    main()
