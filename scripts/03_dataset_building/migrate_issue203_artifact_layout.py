#!/usr/bin/env python3
"""Migrate the two Issue #203 reference artifacts to their canonical layout."""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from sim_swim.analysis.issue203_artifact_layout import main


if __name__ == "__main__":
    main()
