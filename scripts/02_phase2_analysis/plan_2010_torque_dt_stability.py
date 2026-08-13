#!/usr/bin/env python3
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parents[2] / "src"))

from sim_swim.analysis.torque_dt_stability import main

if __name__ == "__main__":
    main()
