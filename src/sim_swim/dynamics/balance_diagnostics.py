"""Same-state component force/torque telemetry (diagnostics, not new QC gates)."""

from __future__ import annotations

import numpy as np

COMPONENTS = (
    "spring",
    "bend",
    "torsion",
    "hook",
    "frame",
    "repulsion",
    "motor",
    "total",
)
DRIVE_COLUMNS = (
    "motor_drive_raw_transverse_ratio_max",
    "motor_drive_transverse_ratio_max",
    "motor_drive_axial_error_ratio_max",
    "motor_drive_correction_force_norm_max_N",
    "motor_drive_correction_relative_norm_max",
    "motor_drive_force_norm_max_N",
)
BALANCE_COLUMNS = (
    [
        f"component_{component}_{quantity}_{suffix}"
        for component in COMPONENTS
        for quantity, unit in (("force", "N"), ("torque", "Nm"))
        for suffix in (f"x_{unit}", f"y_{unit}", f"z_{unit}", f"norm_{unit}")
    ]
    + list(DRIVE_COLUMNS)
    + [
        "force_evaluation_t_s",
        "motor_reaction_support_count_min",
        "motor_reaction_support_count_max",
        "motor_reaction_fallback_used",
        "motor_reaction_solver_success",
    ]
)


def component_balance(
    positions: np.ndarray, components: dict[str, np.ndarray]
) -> dict[str, float]:
    arms = positions - positions.mean(axis=0)
    out = {}
    for name, force in components.items():
        for quantity, unit, vector in (
            ("force", "N", force.sum(axis=0)),
            ("torque", "Nm", np.cross(arms, force).sum(axis=0)),
        ):
            for coordinate, value in zip("xyz", vector, strict=True):
                out[f"component_{name}_{quantity}_{coordinate}_{unit}"] = float(value)
            out[f"component_{name}_{quantity}_norm_{unit}"] = float(
                np.linalg.norm(vector)
            )
    return out
