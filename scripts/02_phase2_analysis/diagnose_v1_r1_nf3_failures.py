#!/usr/bin/env python3
"""Diagnose Phase 2 Issue #158 v1 r1 n=3 long-run failures."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.phase2_158_diagnostics import main


if __name__ == "__main__":
    main()
