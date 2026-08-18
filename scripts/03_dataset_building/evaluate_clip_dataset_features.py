#!/usr/bin/env python3
"""Evaluate reusable 3D source and 2D pixel features for a clip dataset."""

from __future__ import annotations
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from flagella_estimation.phase3.feature_comparison import evaluate

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("dataset_dir", type=Path)
parser.add_argument("output_dir", type=Path, nargs="?")
args = parser.parse_args()
print(
    evaluate(
        args.dataset_dir, args.output_dir or args.dataset_dir / "feature_comparison"
    )
)
