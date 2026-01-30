"""3D可視化（スタブ）。"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

from flagella_sim.sim.core import SimulationState


def save_swim_movie(states: Iterable[SimulationState], out_dir: Path) -> None:
    """3D軌跡の可視化動画/最終フレームを保存する（未実装）。"""

    out_dir.mkdir(parents=True, exist_ok=True)
    note = out_dir / "TODO_render3d.txt"
    note.write_text(
        "3D rendering is not implemented yet. This file is a placeholder.",
        encoding="utf-8",
    )
