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
    torque_ramp_enabled: bool = False,
    torque_ramp_duration_s: float = 0.0,
    enable_switching: bool = False,
    run_tau: float = 0.0,
    tumble_tau: float = 0.0,
    reverse_flagellum: bool = False,
) -> list[np.ndarray]:
    """Reconstruct normalized nominal local-twist segment weights for replay."""
    if segment_count <= 0 or dt_s <= 0:
        raise ValueError("segment_count and dt_s must be positive")
    if profile == "uniform":
        weight = nominal_segment_weights(profile, segment_count)
        return [weight.copy() for _ in times_s]
    if profile not in {
        "diffusive",
        "diffusive_sqrt",
        "diffusive_floor_0p2",
        "diffusive_floor_0p4",
    }:
        raise ValueError(f"unsupported composite torque profile: {profile}")

    orientation = np.zeros(segment_count, dtype=float)
    current_step = 0
    result: list[np.ndarray] = []

    def drive_scale(t_s: float) -> float:
        if not torque_ramp_enabled or torque_ramp_duration_s <= 0:
            return 1.0
        x = min(max(t_s / torque_ramp_duration_s, 0.0), 1.0)
        return x * x * (3.0 - 2.0 * x)

    def torque_sign(t_s: float) -> float:
        if not enable_switching or not reverse_flagellum:
            return 1.0
        cycle = max(run_tau + tumble_tau, 1e-12)
        return -1.0 if (t_s / dt_s) * dt_s % cycle >= run_tau else 1.0

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
            step_t_s = current_step * dt_s
            drive_rate = 2.0 * math.pi * 2.2 * (float(torque_Nm) / 2.0e-20)
            orientation[0] += (
                drive_rate * drive_scale(step_t_s) * torque_sign(step_t_s) * dt_s
            )
            current_step += 1
        weight = np.abs(orientation)
        maximum = float(np.max(weight)) if weight.size else 0.0
        weight = np.ones_like(weight) if maximum <= 1.0e-12 else weight / maximum
        if profile == "diffusive_sqrt":
            weight = np.sqrt(weight)
        elif profile == "diffusive_floor_0p2":
            weight = 0.2 + 0.8 * weight
        elif profile == "diffusive_floor_0p4":
            weight = 0.4 + 0.6 * weight
        result.append(weight / max(float(np.sum(weight)), 1.0e-12))
    return result
