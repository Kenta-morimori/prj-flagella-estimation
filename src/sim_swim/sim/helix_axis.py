"""Flagellar helix-axis diagnostics."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class BodyAxisEstimate:
    """Body axis and rear direction estimate."""

    axis: np.ndarray
    rear_direction: np.ndarray


@dataclass(frozen=True)
class HelixAxisEstimate:
    """Per-flagellum helix-axis estimate."""

    flag_id: int
    origin: np.ndarray
    axis: np.ndarray
    line_start: np.ndarray
    line_end: np.ndarray
    fit_r2: float
    degenerate: bool


def estimate_body_axis(
    positions: np.ndarray,
    body_layer_indices: list[np.ndarray],
    body_indices: np.ndarray | None = None,
) -> BodyAxisEstimate:
    """Estimate body long axis and rear direction from body bead layers."""

    axis = np.zeros(3, dtype=float)
    if len(body_layer_indices) >= 2:
        first = body_layer_indices[0].astype(int, copy=False)
        last = body_layer_indices[-1].astype(int, copy=False)
        c_first = np.mean(positions[first], axis=0)
        c_last = np.mean(positions[last], axis=0)
        axis = c_last - c_first
    elif body_indices is not None and body_indices.size >= 2:
        i0 = int(body_indices[0])
        i1 = int(body_indices[-1])
        axis = positions[i1] - positions[i0]

    norm = float(np.linalg.norm(axis))
    if norm <= 1.0e-18:
        axis = np.array([1.0, 0.0, 0.0], dtype=float)
    else:
        axis = axis / norm
    return BodyAxisEstimate(axis=axis, rear_direction=-axis)


def estimate_flag_helix_axis(
    positions: np.ndarray,
    flag_indices: np.ndarray,
    flag_id: int,
) -> HelixAxisEstimate:
    """Estimate a flagellar helix center axis from beads after the hook bead.

    ``flag_indices[0]`` is treated as the hook-side first bead and is excluded.
    The axis is estimated from ``flag_indices[1:]`` so hook length drift does not
    bias the helix-axis direction.
    """

    idx = np.asarray(flag_indices, dtype=int)
    helix_idx = idx[1:]
    if helix_idx.size < 2:
        nan_vec = np.full(3, float("nan"), dtype=float)
        return HelixAxisEstimate(
            flag_id=int(flag_id),
            origin=nan_vec.copy(),
            axis=nan_vec.copy(),
            line_start=nan_vec.copy(),
            line_end=nan_vec.copy(),
            fit_r2=float("nan"),
            degenerate=True,
        )

    pts = positions[helix_idx]
    origin = np.mean(pts, axis=0)
    centered = pts - origin
    _, singular_values, vh = np.linalg.svd(centered, full_matrices=False)
    variance = singular_values * singular_values
    variance_sum = float(np.sum(variance))
    if variance_sum <= 1.0e-30:
        axis = np.array([1.0, 0.0, 0.0], dtype=float)
        fit_r2 = float("nan")
        degenerate = True
    else:
        axis = vh[0]
        axis = axis / max(float(np.linalg.norm(axis)), 1.0e-18)
        fit_r2 = float(variance[0] / variance_sum)
        degenerate = False

    direction_hint = pts[-1] - pts[0]
    if float(np.dot(axis, direction_hint)) < 0.0:
        axis = -axis

    projections = centered @ axis
    line_start = origin + float(np.min(projections)) * axis
    line_end = origin + float(np.max(projections)) * axis
    return HelixAxisEstimate(
        flag_id=int(flag_id),
        origin=origin,
        axis=axis,
        line_start=line_start,
        line_end=line_end,
        fit_r2=fit_r2,
        degenerate=degenerate,
    )


def angle_deg_between(a: np.ndarray, b: np.ndarray) -> float:
    """Return the unsigned angle between two 3D directions in degrees."""

    norm_a = float(np.linalg.norm(a))
    norm_b = float(np.linalg.norm(b))
    if norm_a <= 1.0e-18 or norm_b <= 1.0e-18:
        return float("nan")
    cos_angle = float(np.clip(np.dot(a, b) / (norm_a * norm_b), -1.0, 1.0))
    return float(np.rad2deg(np.arccos(cos_angle)))
