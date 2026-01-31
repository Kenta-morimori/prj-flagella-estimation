from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class OutputPaths:
    root: Path
    tracking_dir: Path
    sim_dir: Path
    render_dir: Path
    render2d_dir: Path


def make_output_dirs(base_dir: str | Path, date: str, time: str) -> OutputPaths:
    root = Path(base_dir) / date / time
    tracking_dir = root / "tracking"
    sim_dir = root / "sim"
    render_dir = root / "render"
    render2d_dir = root / "render2d"

    # 二重実行の上書きを避けるため exist_ok=False で作成する
    tracking_dir.mkdir(parents=True, exist_ok=False)
    sim_dir.mkdir(parents=True, exist_ok=False)
    render_dir.mkdir(parents=True, exist_ok=False)
    render2d_dir.mkdir(parents=True, exist_ok=False)

    return OutputPaths(
        root=root,
        tracking_dir=tracking_dir,
        sim_dir=sim_dir,
        render_dir=render_dir,
        render2d_dir=render2d_dir,
    )
