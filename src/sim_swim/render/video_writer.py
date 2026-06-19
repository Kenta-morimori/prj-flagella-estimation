"""Shared MP4 writer helpers for render outputs."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import cv2

DEFAULT_MP4_CODECS: tuple[str, ...] = ("avc1", "H264", "mp4v")


@dataclass(frozen=True)
class VideoWriterSelection:
    writer: Any
    selected_codec: str
    attempted_codecs: tuple[str, ...]


@dataclass(frozen=True)
class VideoRenderResult:
    path: str
    selected_codec: str
    attempted_codecs: tuple[str, ...]
    fps: float
    frame_size: tuple[int, int]
    frame_count: int

    def to_manifest(self) -> dict[str, Any]:
        return {
            "path": self.path,
            "selected_codec": self.selected_codec,
            "attempted_codecs": list(self.attempted_codecs),
            "fps": self.fps,
            "frame_size": list(self.frame_size),
            "frame_count": self.frame_count,
        }


def open_mp4_writer(
    path: Path,
    *,
    fps: float,
    frame_size: tuple[int, int],
    codec_candidates: Sequence[str] = DEFAULT_MP4_CODECS,
) -> VideoWriterSelection:
    """Open an MP4 writer, preferring H.264 but falling back to mp4v."""

    attempted_codecs = tuple(codec_candidates)
    for codec in attempted_codecs:
        writer = cv2.VideoWriter(
            str(path),
            cv2.VideoWriter_fourcc(*codec),
            fps,
            frame_size,
        )
        if writer.isOpened():
            return VideoWriterSelection(
                writer=writer,
                selected_codec=codec,
                attempted_codecs=attempted_codecs,
            )
        writer.release()
        try:
            path.unlink(missing_ok=True)
        except OSError:
            pass

    raise RuntimeError(
        "Failed to open MP4 writer with codecs: " + ", ".join(attempted_codecs)
    )
