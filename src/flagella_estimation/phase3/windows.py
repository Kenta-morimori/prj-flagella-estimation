"""Window generation for Phase 3 clip datasets."""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class FrameWindow:
    """Inclusive/exclusive source frame window."""

    start: int
    end: int

    @property
    def frame_count(self) -> int:
        return self.end - self.start


def clip_frame_count(duration_s: float, frame_rate_hz: float) -> int:
    """Return the required frame count for a duration at a frame rate."""

    if duration_s <= 0.0:
        raise ValueError("duration_s must be > 0")
    if frame_rate_hz <= 0.0:
        raise ValueError("frame_rate_hz must be > 0")
    return max(1, int(math.ceil(duration_s * frame_rate_hz)))


def generate_windows(
    *,
    source_frame_count: int,
    frame_rate_hz: float,
    duration_s: float,
    policy: str = "non_overlap",
    overlap_stride_fraction: float = 0.5,
) -> list[FrameWindow]:
    """Generate fixed-length frame windows.

    Partial windows are intentionally dropped so that every clip has a stable
    tensor shape.
    """

    if source_frame_count <= 0:
        return []
    frames_per_clip = clip_frame_count(duration_s, frame_rate_hz)
    if frames_per_clip > source_frame_count:
        return []
    if policy == "full_run":
        return [FrameWindow(0, source_frame_count)]
    if policy == "non_overlap":
        stride = frames_per_clip
    elif policy == "overlap":
        if not 0.0 < overlap_stride_fraction <= 1.0:
            raise ValueError("overlap_stride_fraction must be in (0, 1]")
        stride = max(1, int(round(frames_per_clip * overlap_stride_fraction)))
    else:
        raise ValueError(f"Unsupported window policy: {policy}")

    return [
        FrameWindow(start, start + frames_per_clip)
        for start in range(0, source_frame_count - frames_per_clip + 1, stride)
    ]
