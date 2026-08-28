"""Shared nominal local-segment torque-weight reconstruction for replay."""

from __future__ import annotations

import math

import numpy as np


def nominal_segment_weights(profile: str, segment_count: int) -> np.ndarray:
    """Return normalized nominal weights on bead-to-bead local segments."""
    if segment_count <= 0:
        raise ValueError("segment_count must be positive")
    if profile != "uniform":
        raise ValueError("nominal_segment_weights is only defined for uniform")
    return np.full(segment_count, 1.0 / segment_count, dtype=float)


def reconstructed_segment_weights(
    profile: str,
    segment_count: int,
    *,
    times_s: np.ndarray,
    dt_s: float,
    torque_Nm: float,
) -> list[np.ndarray]:
    """Reconstruct normalized nominal local-twist segment weights for replay."""
    if segment_count <= 0 or dt_s <= 0:
        raise ValueError("segment_count and dt_s must be positive")
    if profile == "uniform":
        weight = nominal_segment_weights(profile, segment_count)
        return [weight.copy() for _ in times_s]
    if profile != "diffusive":
        raise ValueError(f"unsupported composite torque profile: {profile}")

    orientation = np.zeros(segment_count, dtype=float)
    current_step = 0
    result: list[np.ndarray] = []
    drive_rate = 2.0 * math.pi * 2.2 * (float(torque_Nm) / 2.0e-20)
    for t_s in np.asarray(times_s, dtype=float):
        target_step = max(current_step, int(round(float(t_s) / dt_s)))
        while current_step < target_step:
            lap = np.zeros_like(orientation)
            if segment_count > 1:
                lap[0] = orientation[1] - orientation[0]
                lap[-1] = orientation[-2] - orientation[-1]
                if segment_count > 2:
                    lap[1:-1] = (
                        orientation[:-2] - 2.0 * orientation[1:-1] + orientation[2:]
                    )
            orientation += dt_s * (80.0 * lap - 0.05 * orientation)
            orientation[0] += drive_rate * dt_s
            current_step += 1
        weight = np.abs(orientation)
        maximum = float(np.max(weight)) if weight.size else 0.0
        if maximum <= 1.0e-12:
            weight = np.ones_like(weight)
        else:
            weight = weight / maximum
        result.append(weight / max(float(np.sum(weight)), 1.0e-12))
    return result
