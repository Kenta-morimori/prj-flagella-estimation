"""Body-only 2D rendering utilities for pseudo microscopy clips."""

from __future__ import annotations

from dataclasses import dataclass
import math

import cv2
import numpy as np

from sim_swim.sim.core import SimulationState


RENDER_MODE_BODY_CAPSULE_ORTHOGRAPHIC = "body_capsule_orthographic_v1"


@dataclass(frozen=True)
class BodyCapsuleRenderConfig:
    image_size_px: int
    pixel_size_um: float
    body_length_um: float = 2.0
    body_width_um: float = 1.0
    body_intensity: int = 60
    background_intensity: int = 255
    tracking_center: bool = True


@dataclass(frozen=True)
class BodyCapsuleGeometry:
    bbox_xywh_px: tuple[float, float, float, float]
    crop_xywh_px: tuple[float, float, float, float]
    center_xy_px: tuple[float, float]
    body_axis_angle_rad: float | None
    body_length_px: float
    body_width_px: float


def body_axis_3d_from_state(state: SimulationState) -> np.ndarray:
    """Return the rigid body long-axis direction in world coordinates."""

    q = np.asarray(state.quaternion, dtype=float)
    if q.shape == (4,) and np.isfinite(q).all():
        norm = float(np.linalg.norm(q))
        if norm > 1.0e-12:
            x, y, z, w = q / norm
            axis = np.array(
                [
                    1.0 - 2.0 * (y * y + z * z),
                    2.0 * (x * y + z * w),
                    2.0 * (x * z - y * w),
                ],
                dtype=float,
            )
            axis_norm = float(np.linalg.norm(axis))
            if axis_norm > 1.0e-12:
                return axis / axis_norm

    beads = np.asarray(state.bead_positions_um, dtype=float)
    if beads.ndim == 2 and beads.shape[0] >= 2 and beads.shape[1] >= 3:
        centered = beads[:, :3] - np.mean(beads[:, :3], axis=0)
        try:
            _, singular_values, vh = np.linalg.svd(centered, full_matrices=False)
        except np.linalg.LinAlgError:
            singular_values = np.zeros(1, dtype=float)
            vh = np.asarray([[1.0, 0.0, 0.0]], dtype=float)
        if float(singular_values[0]) > 1.0e-12:
            axis = np.asarray(vh[0], dtype=float)
            axis_norm = float(np.linalg.norm(axis))
            if axis_norm > 1.0e-12:
                return axis / axis_norm

    return np.asarray([1.0, 0.0, 0.0], dtype=float)


def body_axis_xy_from_state(state: SimulationState) -> np.ndarray:
    """Return the normalized image-plane direction of the body long axis."""

    axis_xy = body_axis_3d_from_state(state)[:2]
    axis_xy_norm = float(np.linalg.norm(axis_xy))
    if axis_xy_norm <= 1.0e-12:
        return np.asarray([1.0, 0.0], dtype=float)
    return axis_xy / axis_xy_norm


def render_body_capsule_frame(
    state: SimulationState,
    cfg: BodyCapsuleRenderConfig,
) -> tuple[np.ndarray, BodyCapsuleGeometry]:
    """Render one rigid body capsule frame without flagella or body deformation."""

    if cfg.image_size_px <= 0:
        raise ValueError("image_size_px must be > 0")
    if cfg.pixel_size_um <= 0.0:
        raise ValueError("pixel_size_um must be > 0")
    if cfg.body_length_um <= 0.0:
        raise ValueError("body_length_um must be > 0")
    if cfg.body_width_um <= 0.0:
        raise ValueError("body_width_um must be > 0")
    if cfg.body_length_um < cfg.body_width_um:
        raise ValueError("body_length_um must be >= body_width_um")

    image_size = int(cfg.image_size_px)
    frame = np.full(
        (image_size, image_size),
        int(np.clip(cfg.background_intensity, 0, 255)),
        dtype=np.uint8,
    )
    px_per_um = 1.0 / cfg.pixel_size_um
    intrinsic_body_length_px = float(cfg.body_length_um * px_per_um)
    body_width_px = float(cfg.body_width_um * px_per_um)
    center_um = (
        np.asarray(state.position_um[:2], dtype=float)
        if cfg.tracking_center
        else np.zeros(2, dtype=float)
    )
    state_center_um = np.asarray(state.position_um[:2], dtype=float)
    center_px = (state_center_um - center_um) * px_per_um + image_size / 2.0

    axis_xy = body_axis_3d_from_state(state)[:2]
    projection_scale = float(np.clip(np.linalg.norm(axis_xy), 0.0, 1.0))
    projected_centerline_px = (
        intrinsic_body_length_px - body_width_px
    ) * projection_scale
    body_length_px = body_width_px + projected_centerline_px
    angle_is_observable = projected_centerline_px >= 0.5
    if projection_scale > 1.0e-12:
        axis_xy = axis_xy / projection_scale
    else:
        axis_xy = np.asarray([1.0, 0.0], dtype=float)
    angle_rad = (
        float(math.atan2(axis_xy[1], axis_xy[0])) if angle_is_observable else None
    )
    half_axis = axis_xy * (projected_centerline_px / 2.0)
    p0 = center_px - half_axis
    p1 = center_px + half_axis
    thickness = max(1, int(round(body_width_px)))
    color = int(np.clip(cfg.body_intensity, 0, 255))
    if angle_is_observable:
        cv2.line(
            frame,
            tuple(np.rint(p0).astype(int)),
            tuple(np.rint(p1).astype(int)),
            color,
            thickness,
            cv2.LINE_AA,
        )
    else:
        cv2.circle(
            frame,
            tuple(np.rint(center_px).astype(int)),
            max(1, int(round(body_width_px / 2.0))),
            color,
            -1,
            cv2.LINE_AA,
        )

    radius = body_width_px / 2.0
    min_xy = np.minimum(p0, p1) - radius
    max_xy = np.maximum(p0, p1) + radius
    bbox_xywh = (
        float(min_xy[0]),
        float(min_xy[1]),
        float(max(max_xy[0] - min_xy[0], 1.0)),
        float(max(max_xy[1] - min_xy[1], 1.0)),
    )
    geometry = BodyCapsuleGeometry(
        bbox_xywh_px=bbox_xywh,
        crop_xywh_px=(0.0, 0.0, float(image_size), float(image_size)),
        center_xy_px=(float(center_px[0]), float(center_px[1])),
        body_axis_angle_rad=angle_rad,
        body_length_px=body_length_px,
        body_width_px=body_width_px,
    )
    return frame, geometry


def render_body_capsule_clip(
    states: list[SimulationState],
    cfg: BodyCapsuleRenderConfig,
) -> tuple[np.ndarray, list[BodyCapsuleGeometry]]:
    """Render states into a grayscale uint8 ``(T, H, W)`` body-only clip."""

    if not states:
        raise ValueError("Cannot render an empty clip")
    frames: list[np.ndarray] = []
    geometries: list[BodyCapsuleGeometry] = []
    for state in states:
        frame, geometry = render_body_capsule_frame(state, cfg)
        frames.append(frame)
        geometries.append(geometry)
    return np.stack(frames, axis=0), geometries
