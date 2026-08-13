#!/usr/bin/env python3
"""Issue #183 reference torque campaign plan CLI (simulation は起動しない)。"""

from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[2] / "src"))

from sim_swim.analysis.reference_torque_comparison import main


if __name__ == "__main__":
    main()
