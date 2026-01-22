from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class OutputPaths:
    root: Path
    tracking_dir: Path


def make_output_dirs(base_dir: str | Path, date: str, time: str) -> OutputPaths:
    root = Path(base_dir) / date / time
    tracking_dir = root / "tracking"
    tracking_dir.mkdir(parents=True, exist_ok=False)  # 二重実行の上書き防止
    return OutputPaths(root=root, tracking_dir=tracking_dir)
