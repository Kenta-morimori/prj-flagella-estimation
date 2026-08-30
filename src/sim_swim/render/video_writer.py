"""Shared, PowerPoint-compatible H.264 MP4 writer helpers.

All rendered MP4 artifacts use FFmpeg's ``libx264`` encoder. OpenCV's MP4
backends vary by host and can silently fall back to ``mp4v``; that format is
not an acceptable delivery format for this project.
"""

from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import shutil
import subprocess
from typing import Any, Sequence

import numpy as np


DEFAULT_MP4_CODECS: tuple[str, ...] = ("ffmpeg:libx264",)


def resolve_ffmpeg() -> str:
    """Return a usable FFmpeg executable with the required H.264 encoder.

    cs10 does not expose user-local binaries in its default ``PATH``. Check
    its documented installation location before the system path so rendering
    works both in interactive shells and parallel workers.
    """
    candidates = [Path.home() / ".local" / "bin" / "ffmpeg"]
    path_ffmpeg = shutil.which("ffmpeg")
    if path_ffmpeg:
        candidates.append(Path(path_ffmpeg))
    for candidate in candidates:
        if not candidate.is_file() or not os.access(candidate, os.X_OK):
            continue
        completed = subprocess.run(
            [str(candidate), "-hide_banner", "-encoders"],
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode == 0 and "libx264" in completed.stdout:
            return str(candidate)
    raise RuntimeError(
        "H.264 MP4 rendering requires FFmpeg with libx264. "
        "On cs10 run scripts/cs10/setup_environment.sh, then retry."
    )


class _FFmpegVideoWriter:
    """Small ``cv2.VideoWriter``-compatible BGR rawvideo pipe."""

    def __init__(self, path: Path, *, fps: float, frame_size: tuple[int, int]) -> None:
        self.path = path
        self.frame_size = frame_size
        self._released = False
        ffmpeg = resolve_ffmpeg()
        path.parent.mkdir(parents=True, exist_ok=True)
        self._process = subprocess.Popen(
            [
                ffmpeg,
                "-y",
                "-loglevel",
                "error",
                "-f",
                "rawvideo",
                "-pix_fmt",
                "bgr24",
                "-video_size",
                f"{frame_size[0]}x{frame_size[1]}",
                "-framerate",
                str(fps),
                "-i",
                "-",
                "-an",
                # yuv420p requires even image dimensions. Matplotlib grid
                # canvases may be odd-sized (for example 1439 px wide), so
                # pad at encode time instead of failing after the first frame.
                "-vf",
                "pad=ceil(iw/2)*2:ceil(ih/2)*2",
                "-c:v",
                "libx264",
                "-profile:v",
                "high",
                "-vf",
                "pad=ceil(iw/2)*2:ceil(ih/2)*2:color=black",
                "-pix_fmt",
                "yuv420p",
                "-movflags",
                "+faststart",
                str(path),
            ],
            stdin=subprocess.PIPE,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )

    def isOpened(self) -> bool:
        return not self._released and self._process.poll() is None

    def write(self, frame: np.ndarray) -> None:
        expected = (self.frame_size[1], self.frame_size[0], 3)
        if frame.shape != expected or frame.dtype != np.uint8:
            raise ValueError(
                f"MP4 frame must be uint8 BGR with shape {expected}; "
                f"received {frame.shape} {frame.dtype}"
            )
        if not self.isOpened() or self._process.stdin is None:
            stderr = ""
            if self._process.stderr is not None:
                stderr = self._process.stderr.read().decode("utf-8", errors="replace")
            raise RuntimeError(
                "FFmpeg MP4 writer is not open "
                f"(exit_code={self._process.poll()}): {stderr.strip()}"
            )
        self._process.stdin.write(np.ascontiguousarray(frame).tobytes())

    def release(self) -> None:
        if self._released:
            return
        self._released = True
        if self._process.stdin is not None:
            self._process.stdin.close()
        exit_code = self._process.wait()
        if exit_code != 0:
            stderr = ""
            if self._process.stderr is not None:
                stderr = self._process.stderr.read().decode("utf-8", errors="replace")
            self.path.unlink(missing_ok=True)
            raise RuntimeError(
                "FFmpeg failed while encoding H.264 MP4 "
                f"(exit_code={exit_code}): {self.path}: {stderr.strip()}"
            )


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
    """Open a writer that always emits H.264/high-profile/yuv420p MP4.

    ``codec_candidates`` remains as a compatibility parameter for callers but
    no longer authorizes fallback to a non-H.264 codec.
    """
    attempted_codecs = tuple(codec_candidates)
    if attempted_codecs != DEFAULT_MP4_CODECS:
        raise ValueError(
            "Only H.264 output is supported; codec_candidates must be "
            f"{DEFAULT_MP4_CODECS}"
        )
    if fps <= 0 or frame_size[0] <= 0 or frame_size[1] <= 0:
        raise ValueError("fps and frame_size must be positive")
    writer = _FFmpegVideoWriter(path, fps=fps, frame_size=frame_size)
    if not writer.isOpened():
        writer.release()
        raise RuntimeError(f"Failed to open FFmpeg H.264 writer: {path}")
    return VideoWriterSelection(
        writer=writer,
        selected_codec="libx264",
        attempted_codecs=DEFAULT_MP4_CODECS,
    )
