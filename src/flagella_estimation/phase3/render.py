"""Lightweight Phase 2 state archive to Phase 3 clip rendering."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from sim_swim.sim.core import SimulationState


@dataclass(frozen=True)
class FrameGeometry:
    bbox_xywh_px: tuple[float, float, float, float]
    crop_xywh_px: tuple[float, float, float, float]
    center_xy_px: tuple[float, float]
    body_axis_angle_rad: float | None
    body_length_px: float | None
    body_width_px: float | None


def select_frames(
    states: list[SimulationState], frame_rate_hz: float
) -> list[SimulationState]:
    """Select approximately evenly-spaced output frames from raw simulation states."""

    if not states:
        return []
    if frame_rate_hz <= 0.0:
        raise ValueError("frame_rate_hz must be > 0")
    interval_s = 1.0 / frame_rate_hz
    selected: list[SimulationState] = []
    next_t = float(states[0].t)
    for state in states:
        if float(state.t) + 1.0e-12 >= next_t:
            selected.append(state)
            next_t += interval_s
    if selected and selected[-1] is not states[-1]:
        selected.append(states[-1])
    return selected


def render_state_frame(
    state: SimulationState,
    *,
    image_size_px: int,
    pixel_size_um: float,
) -> tuple[np.ndarray, FrameGeometry]:
    """Render archived bead positions as a centered grayscale numpy frame."""

    if image_size_px <= 0:
        raise ValueError("image_size_px must be > 0")
    if pixel_size_um <= 0.0:
        raise ValueError("pixel_size_um must be > 0")

    frame = np.full((image_size_px, image_size_px), 255, dtype=np.uint8)
    beads = np.asarray(state.bead_positions_um, dtype=float)
    if beads.size == 0:
        center = (image_size_px / 2.0, image_size_px / 2.0)
        geom = FrameGeometry(
            bbox_xywh_px=(center[0], center[1], 0.0, 0.0),
            crop_xywh_px=(0.0, 0.0, float(image_size_px), float(image_size_px)),
            center_xy_px=center,
            body_axis_angle_rad=None,
            body_length_px=None,
            body_width_px=None,
        )
        return frame, geom

    xy_um = beads[:, :2]
    center_um = np.asarray(state.position_um[:2], dtype=float)
    xy_px = (xy_um - center_um) / pixel_size_um + image_size_px / 2.0
    rounded = np.rint(xy_px).astype(int)
    valid = (
        (rounded[:, 0] >= 0)
        & (rounded[:, 0] < image_size_px)
        & (rounded[:, 1] >= 0)
        & (rounded[:, 1] < image_size_px)
    )
    for x_px, y_px in rounded[valid]:
        x0 = max(int(x_px) - 1, 0)
        x1 = min(int(x_px) + 2, image_size_px)
        y0 = max(int(y_px) - 1, 0)
        y1 = min(int(y_px) + 2, image_size_px)
        frame[y0:y1, x0:x1] = 60

    min_xy = np.min(xy_px, axis=0)
    max_xy = np.max(xy_px, axis=0)
    bbox_w, bbox_h = np.maximum(max_xy - min_xy, 1.0)
    center_xy = np.mean(xy_px, axis=0)
    geom = FrameGeometry(
        bbox_xywh_px=(float(min_xy[0]), float(min_xy[1]), float(bbox_w), float(bbox_h)),
        crop_xywh_px=(0.0, 0.0, float(image_size_px), float(image_size_px)),
        center_xy_px=(float(center_xy[0]), float(center_xy[1])),
        body_axis_angle_rad=None,
        body_length_px=float(max(bbox_w, bbox_h)),
        body_width_px=float(min(bbox_w, bbox_h)),
    )
    return frame, geom


def render_clip_array(
    states: list[SimulationState],
    *,
    image_size_px: int,
    pixel_size_um: float,
) -> tuple[np.ndarray, list[FrameGeometry]]:
    """Render a list of states into a uint8 `(T, H, W)` clip array."""

    frames: list[np.ndarray] = []
    geometries: list[FrameGeometry] = []
    for state in states:
        frame, geometry = render_state_frame(
            state,
            image_size_px=image_size_px,
            pixel_size_um=pixel_size_um,
        )
        frames.append(frame)
        geometries.append(geometry)
    if not frames:
        raise ValueError("Cannot render an empty clip")
    return np.stack(frames, axis=0), geometries
