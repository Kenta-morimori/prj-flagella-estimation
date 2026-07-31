"""Lightweight Phase 2 state archive to Phase 3 clip rendering."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from sim_swim.render.body2d import (
    RENDER_MODE_BODY_CAPSULE_ORTHOGRAPHIC,
    BodyCapsuleRenderConfig,
    render_body_capsule_clip,
    render_body_capsule_frame,
)
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
    body_length_um: float = 2.0,
    body_width_um: float = 1.0,
    body_intensity: int = 60,
) -> tuple[np.ndarray, FrameGeometry]:
    """Render a centered rigid body capsule grayscale frame."""

    frame, capsule_geom = render_body_capsule_frame(
        state,
        BodyCapsuleRenderConfig(
            image_size_px=image_size_px,
            pixel_size_um=pixel_size_um,
            body_length_um=body_length_um,
            body_width_um=body_width_um,
            body_intensity=body_intensity,
            tracking_center=True,
        ),
    )
    geom = FrameGeometry(
        bbox_xywh_px=capsule_geom.bbox_xywh_px,
        crop_xywh_px=capsule_geom.crop_xywh_px,
        center_xy_px=capsule_geom.center_xy_px,
        body_axis_angle_rad=capsule_geom.body_axis_angle_rad,
        body_length_px=capsule_geom.body_length_px,
        body_width_px=capsule_geom.body_width_px,
    )
    return frame, geom


def render_clip_array(
    states: list[SimulationState],
    *,
    image_size_px: int,
    pixel_size_um: float,
    body_length_um: float = 2.0,
    body_width_um: float = 1.0,
    body_intensity: int = 60,
) -> tuple[np.ndarray, list[FrameGeometry]]:
    """Render a list of states into a uint8 `(T, H, W)` clip array."""

    clip, capsule_geometries = render_body_capsule_clip(
        states,
        BodyCapsuleRenderConfig(
            image_size_px=image_size_px,
            pixel_size_um=pixel_size_um,
            body_length_um=body_length_um,
            body_width_um=body_width_um,
            body_intensity=body_intensity,
            tracking_center=True,
        ),
    )
    geometries = [
        FrameGeometry(
            bbox_xywh_px=geometry.bbox_xywh_px,
            crop_xywh_px=geometry.crop_xywh_px,
            center_xy_px=geometry.center_xy_px,
            body_axis_angle_rad=geometry.body_axis_angle_rad,
            body_length_px=geometry.body_length_px,
            body_width_px=geometry.body_width_px,
        )
        for geometry in capsule_geometries
    ]
    return clip, geometries


RENDER_MODE = RENDER_MODE_BODY_CAPSULE_ORTHOGRAPHIC
