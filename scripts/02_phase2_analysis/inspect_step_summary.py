#!/usr/bin/env python3
"""run_summary.json を入口に、step_summary.csv の限定区間だけを表示する。"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.step_summary_inspection import main


if __name__ == "__main__":
    main()
