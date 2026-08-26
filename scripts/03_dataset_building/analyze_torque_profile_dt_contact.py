#!/usr/bin/env python3
"""Analyze reusable #203 references plus a compact n=3 dt diagnostic campaign."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from sim_swim.analysis.torque_profile_dt_contact import main


if __name__ == "__main__":
    main()
