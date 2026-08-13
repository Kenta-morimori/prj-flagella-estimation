#!/usr/bin/env python3
"""Analyze Issue #163 spring-formulation comparison results."""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.spring_formulations_workflow import main


if __name__ == "__main__":
    main()
