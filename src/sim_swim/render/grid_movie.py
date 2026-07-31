"""Shared helpers for MP4 grid replay artifacts."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import cv2
import numpy as np

from sim_swim.render.video_writer import VideoRenderResult, open_mp4_writer


@dataclass(frozen=True)
class GridLayout:
    rows: int
    cols: int
    cell_width_px: int
    cell_height_px: int

    @property
    def frame_size(self) -> tuple[int, int]:
        return (self.cell_width_px * self.cols, self.cell_height_px * self.rows)


def auto_grid_layout(
    item_count: int,
    *,
    cell_width_px: int,
    cell_height_px: int,
    max_cols: int = 2,
) -> GridLayout:
    if item_count <= 0:
        raise ValueError("item_count must be positive")
    cols = min(max(1, int(max_cols)), item_count)
    rows = int(np.ceil(item_count / cols))
    return GridLayout(
        rows=rows,
        cols=cols,
        cell_width_px=int(cell_width_px),
        cell_height_px=int(cell_height_px),
    )


def compose_grid_frame(
    panels: list[np.ndarray],
    layout: GridLayout,
    *,
    background_bgr: tuple[int, int, int] = (255, 255, 255),
) -> np.ndarray:
    frame_w, frame_h = layout.frame_size
    canvas = np.full((frame_h, frame_w, 3), background_bgr, dtype=np.uint8)
    for panel_index, panel in enumerate(panels):
        row = panel_index // layout.cols
        col = panel_index % layout.cols
        if row >= layout.rows:
            break
        y0 = row * layout.cell_height_px
        x0 = col * layout.cell_width_px
        if panel.shape[:2] != (layout.cell_height_px, layout.cell_width_px):
            panel = cv2.resize(
                panel,
                (layout.cell_width_px, layout.cell_height_px),
                interpolation=cv2.INTER_AREA,
            )
        canvas[y0 : y0 + layout.cell_height_px, x0 : x0 + layout.cell_width_px] = panel
    return canvas


def write_mp4_grid(
    path: Path,
    *,
    frames: Iterable[np.ndarray],
    fps: float,
) -> VideoRenderResult:
    writer_selection = None
    writer = None
    frame_count = 0
    last_frame: np.ndarray | None = None
    for frame in frames:
        if writer is None:
            writer_selection = open_mp4_writer(
                path,
                fps=max(float(fps), 1.0e-6),
                frame_size=(frame.shape[1], frame.shape[0]),
            )
            writer = writer_selection.writer
        writer.write(frame)
        frame_count += 1
        last_frame = frame
    if writer is not None:
        writer.release()
    if writer_selection is None or last_frame is None:
        raise RuntimeError("No frames were provided to write_mp4_grid")
    return VideoRenderResult(
        path=str(path),
        selected_codec=writer_selection.selected_codec,
        attempted_codecs=writer_selection.attempted_codecs,
        fps=float(fps),
        frame_size=(last_frame.shape[1], last_frame.shape[0]),
        frame_count=frame_count,
    )
