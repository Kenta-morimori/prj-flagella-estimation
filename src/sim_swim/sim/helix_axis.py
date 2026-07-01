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


@dataclass(frozen=True)
class HelixAxisAlignmentMetrics:
    """Aggregate alignment metrics across multiple flagellar helix axes."""

    pair_angle_deg_mean: float
    pair_angle_deg_max: float
    mean_deviation_deg_max: float
    alignment_order: float


@dataclass(frozen=True)
class HelixAxisCenteredMetrics:
    """Metrics for helix shape around its estimated center axis."""

    phase_deg: float
    fit_r2: float
    radius_mean_m: float
    radius_std_m: float
    radius_cv: float
    root_offset_m: float
    degenerate: bool


def _degenerate_helix_axis_estimate(flag_id: int) -> HelixAxisEstimate:
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


def _degenerate_centered_metrics() -> HelixAxisCenteredMetrics:
    return HelixAxisCenteredMetrics(
        phase_deg=float("nan"),
        fit_r2=float("nan"),
        radius_mean_m=float("nan"),
        radius_std_m=float("nan"),
        radius_cv=float("nan"),
        root_offset_m=float("nan"),
        degenerate=True,
    )


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
        return _degenerate_helix_axis_estimate(flag_id)

    pts = positions[helix_idx]
    if not np.isfinite(pts).all():
        return _degenerate_helix_axis_estimate(flag_id)

    origin = np.mean(pts, axis=0)
    centered = pts - origin
    try:
        _, singular_values, vh = np.linalg.svd(centered, full_matrices=False)
    except np.linalg.LinAlgError:
        return _degenerate_helix_axis_estimate(flag_id)

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


def helix_axis_centered_metrics(
    positions: np.ndarray,
    flag_indices: np.ndarray,
    estimate: HelixAxisEstimate,
    reference_direction: np.ndarray,
) -> HelixAxisCenteredMetrics:
    """Measure how cleanly helix beads wrap around the estimated center axis."""

    idx = np.asarray(flag_indices, dtype=int)
    helix_idx = idx[1:]
    if helix_idx.size < 3 or estimate.degenerate:
        return _degenerate_centered_metrics()

    pts = positions[helix_idx]
    axis = np.asarray(estimate.axis, dtype=float)
    origin = np.asarray(estimate.origin, dtype=float)
    ref = np.asarray(reference_direction, dtype=float)
    if (
        not np.isfinite(pts).all()
        or not np.isfinite(axis).all()
        or not np.isfinite(origin).all()
        or not np.isfinite(ref).all()
    ):
        return _degenerate_centered_metrics()

    axis_norm = float(np.linalg.norm(axis))
    if axis_norm <= 1.0e-18:
        return _degenerate_centered_metrics()
    axis = axis / axis_norm

    e1 = ref - float(np.dot(ref, axis)) * axis
    if float(np.linalg.norm(e1)) <= 1.0e-18:
        rel0 = pts[0] - origin
        e1 = rel0 - float(np.dot(rel0, axis)) * axis
    if float(np.linalg.norm(e1)) <= 1.0e-18:
        fallback = np.array([1.0, 0.0, 0.0], dtype=float)
        if abs(float(np.dot(fallback, axis))) > 0.9:
            fallback = np.array([0.0, 1.0, 0.0], dtype=float)
        e1 = fallback - float(np.dot(fallback, axis)) * axis
    e1_norm = float(np.linalg.norm(e1))
    if e1_norm <= 1.0e-18:
        return _degenerate_centered_metrics()
    e1 = e1 / e1_norm
    e2 = np.cross(axis, e1)
    e2 = e2 / max(float(np.linalg.norm(e2)), 1.0e-18)

    rel = pts - origin
    axial = rel @ axis
    radial = rel - np.outer(axial, axis)
    radius = np.linalg.norm(radial, axis=1)
    radius_mean = float(np.mean(radius))
    radius_std = float(np.std(radius))
    radius_cv = radius_std / max(radius_mean, 1.0e-18)

    u = radial @ e1
    v = radial @ e2
    theta = np.unwrap(np.arctan2(v, u))
    if theta.size < 3 or not np.isfinite(theta).all():
        return _degenerate_centered_metrics()
    slope, intercept = np.polyfit(axial, theta, 1)
    pred = slope * axial + intercept
    ss_res = float(np.sum((theta - pred) ** 2))
    ss_tot = float(np.sum((theta - float(np.mean(theta))) ** 2))
    fit_r2 = 1.0 - ss_res / max(ss_tot, 1.0e-30)

    root = positions[int(idx[0])]
    root_rel = root - origin
    root_radial = root_rel - float(np.dot(root_rel, axis)) * axis
    root_offset = float(np.linalg.norm(root_radial))

    return HelixAxisCenteredMetrics(
        phase_deg=float(np.rad2deg(intercept)),
        fit_r2=float(fit_r2),
        radius_mean_m=radius_mean,
        radius_std_m=radius_std,
        radius_cv=float(radius_cv),
        root_offset_m=root_offset,
        degenerate=False,
    )


def helix_axis_alignment_metrics(axes: list[np.ndarray]) -> HelixAxisAlignmentMetrics:
    """Return pairwise and mean-axis alignment metrics for oriented axes."""

    unit_axes: list[np.ndarray] = []
    for axis_raw in axes:
        axis = np.asarray(axis_raw, dtype=float)
        norm = float(np.linalg.norm(axis))
        if norm <= 1.0e-18 or not np.isfinite(axis).all():
            continue
        unit_axes.append(axis / norm)

    if not unit_axes:
        return HelixAxisAlignmentMetrics(
            pair_angle_deg_mean=float("nan"),
            pair_angle_deg_max=float("nan"),
            mean_deviation_deg_max=float("nan"),
            alignment_order=float("nan"),
        )
    if len(unit_axes) == 1:
        return HelixAxisAlignmentMetrics(
            pair_angle_deg_mean=0.0,
            pair_angle_deg_max=0.0,
            mean_deviation_deg_max=0.0,
            alignment_order=1.0,
        )

    pair_angles: list[float] = []
    for i, axis_i in enumerate(unit_axes):
        for axis_j in unit_axes[i + 1 :]:
            pair_angles.append(angle_deg_between(axis_i, axis_j))

    mean_vec = np.mean(np.asarray(unit_axes, dtype=float), axis=0)
    alignment_order = float(np.linalg.norm(mean_vec))
    mean_norm = float(np.linalg.norm(mean_vec))
    if mean_norm <= 1.0e-18:
        mean_deviation_max = 180.0
    else:
        mean_axis = mean_vec / mean_norm
        mean_deviation_max = float(
            np.max([angle_deg_between(axis, mean_axis) for axis in unit_axes])
        )

    pair_arr = np.asarray(pair_angles, dtype=float)
    return HelixAxisAlignmentMetrics(
        pair_angle_deg_mean=float(np.nanmean(pair_arr)),
        pair_angle_deg_max=float(np.nanmax(pair_arr)),
        mean_deviation_deg_max=mean_deviation_max,
        alignment_order=alignment_order,
    )
