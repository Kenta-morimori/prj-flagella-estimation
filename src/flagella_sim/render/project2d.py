"""3D軌跡から2D投影フレームを生成するためのスタブ。"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

from flagella_sim.sim.core import SimulationState
from flagella_sim.sim.params import RenderParams


def project_states(
    states: Iterable[SimulationState],
    render_cfg: RenderParams,
    out_dir: Path,
) -> None:
    """軌跡を2Dへ投影しPNG連番・mp4を出力する（現状はスタブ）。

    将来的には orthographic 投影＋描画を実装する。現時点ではディレクトリのみ作成する。
    """

    frames_dir = out_dir / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)
    # 未実装の印としてメモを残す
    note = out_dir / "TODO_render2d.txt"
    note.write_text(
        "render2d generation is not implemented yet. This file is a placeholder.",
        encoding="utf-8",
    )
