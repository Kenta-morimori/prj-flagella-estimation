"""single flagellum の multi-step helix retention gate。"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np

HELIX_RETENTION_BOND_REL_ERR_MAX_LIMIT = 0.25
HELIX_RETENTION_BEND_ERR_MAX_DEG_LIMIT = 30.0
HELIX_RETENTION_TORSION_ERR_MAX_DEG_LIMIT = 60.0
HELIX_RETENTION_HOOK_REL_ERR_MAX_LIMIT = 0.5

REQUIRED_HELIX_RETENTION_NUMERIC_COLUMNS = (
    "flag_phase_deg",
    "flag_phase_rate_hz",
    "flag_helix_spin_phase_deg",
    "flag_helix_spin_rate_hz",
    "flag_helix_spin_fit_r2",
    "hook_len_rel_err_max",
    "local_attach_first_rel_err",
    "flag_bond_rel_err_max",
    "flag_bend_err_max_deg",
    "flag_torsion_err_max_deg",
)

REQUIRED_HELIX_RETENTION_ZERO_COUNTER_COLUMNS = (
    "motor_degenerate_axis_count",
    "motor_bond_length_clipped_count",
)


def _parse_float(value: str | float | None) -> float:
    if value is None:
        return float("nan")
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _parse_bool(value: str | bool | None) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def _empty_summary(category: str, reason: str) -> dict[str, bool | int | float | str]:
    return {
        "helix_retention_pass": False,
        "step_count": 0,
        "first_fail_step": -1,
        "first_fail_category": category,
        "first_fail_reason": reason,
        "median_abs_flag_phase_rate_hz": float("nan"),
        "net_abs_flag_root_revolutions": float("nan"),
        "signed_flag_root_rate_hz": float("nan"),
        "median_abs_flag_helix_spin_rate_hz": float("nan"),
        "net_abs_flag_helix_spin_revolutions": float("nan"),
        "signed_flag_helix_spin_rate_hz": float("nan"),
        "flag_helix_spin_direction_consistency": float("nan"),
        "helix_to_root_net_rotation_ratio": float("nan"),
        "min_flag_helix_spin_fit_r2": float("nan"),
        "max_hook_len_rel_err": float("nan"),
        "max_local_attach_first_rel_err": float("nan"),
        "max_flag_bond_rel_err": float("nan"),
        "max_flag_bend_err_deg": float("nan"),
        "max_flag_torsion_err_deg": float("nan"),
    }


def summarize_single_flagellum_helix_retention(
    rows: Sequence[Mapping[str, str | float | None]],
    *,
    skip_initial_steps: int = 1,
    min_steps: int = 50,
    min_median_abs_spin_rate_hz: float = 1.0,
    min_net_abs_spin_revolutions: float = 1.0,
    min_direction_consistency: float = 0.5,
    min_helix_spin_fit_r2: float = 0.5,
) -> dict[str, bool | int | float | str]:
    """step summary rows から multi-step helix retention gate を判定する。

    Phase 2.6 では body diagnostics が長時間では出力されないため、この gate は
    `step_summary.csv` の non-body / motor / flagellum 診断だけを正本とする。
    """

    if not rows:
        return _empty_summary("empty", "step_summary rows are empty")

    eval_rows = list(rows[int(skip_initial_steps) :])
    if len(eval_rows) < min_steps:
        summary = _empty_summary(
            "insufficient_steps",
            f"step_count={len(eval_rows)} < min_steps={min_steps}",
        )
        summary["step_count"] = len(eval_rows)
        return summary

    first_fail_step = -1
    first_fail_category = "none"
    first_fail_reason = "none"
    for row_idx, row in enumerate(eval_rows):
        step_value = _parse_float(row.get("step"))
        step = int(step_value) if np.isfinite(step_value) else row_idx

        if not _parse_bool(row.get("finite_pass")):
            first_fail_step = step
            first_fail_category = "finite"
            first_fail_reason = "finite_pass=False"
            break
        if not _parse_bool(row.get("shape_pass_nonbody")):
            first_fail_step = step
            first_fail_category = str(row.get("first_fail_category_nonbody", "nonbody"))
            first_fail_reason = "shape_pass_nonbody=False"
            break

        for col in REQUIRED_HELIX_RETENTION_NUMERIC_COLUMNS:
            value = _parse_float(row.get(col))
            if not np.isfinite(value):
                first_fail_step = step
                first_fail_category = "nonfinite"
                first_fail_reason = f"non-finite value: {col}"
                break
        if first_fail_category != "none":
            break

        for col in REQUIRED_HELIX_RETENTION_ZERO_COUNTER_COLUMNS:
            count_value = _parse_float(row.get(col))
            if not np.isfinite(count_value):
                first_fail_step = step
                first_fail_category = "nonfinite"
                first_fail_reason = f"non-finite counter: {col}"
                break
            count = int(count_value)
            if count > 0:
                first_fail_step = step
                first_fail_category = "motor_split"
                first_fail_reason = f"counter incremented: {col}={count}"
                break
        if first_fail_category != "none":
            break

    def series(name: str) -> np.ndarray:
        return np.asarray([_parse_float(row.get(name)) for row in eval_rows])

    def max_col(name: str) -> float:
        values = series(name)
        return float(np.max(values)) if np.isfinite(values).all() else float("nan")

    max_hook = max_col("hook_len_rel_err_max")
    max_attach = max_col("local_attach_first_rel_err")
    max_bond = max_col("flag_bond_rel_err_max")
    max_bend = max_col("flag_bend_err_max_deg")
    max_torsion = max_col("flag_torsion_err_max_deg")
    root_phases = series("flag_phase_deg")
    phase_rates = series("flag_phase_rate_hz")
    spin_phases = series("flag_helix_spin_phase_deg")
    spin_rates = series("flag_helix_spin_rate_hz")
    spin_fit_r2 = series("flag_helix_spin_fit_r2")
    median_abs_phase_rate_hz = (
        float(np.median(np.abs(phase_rates)))
        if np.isfinite(phase_rates).all()
        else float("nan")
    )
    net_root_deg = float(root_phases[-1] - root_phases[0])
    net_abs_root_revolutions = abs(net_root_deg) / 360.0
    duration_s = float(_parse_float(eval_rows[-1].get("t_s"))) - float(
        _parse_float(eval_rows[0].get("t_s"))
    )
    signed_root_rate_hz = (
        net_root_deg / 360.0 / duration_s
        if np.isfinite(duration_s) and duration_s > 0.0
        else float("nan")
    )
    median_abs_spin_rate_hz = float(np.median(np.abs(spin_rates)))
    spin_phase_deltas = np.diff(spin_phases)
    net_spin_deg = float(spin_phases[-1] - spin_phases[0])
    total_abs_spin_deg = float(np.sum(np.abs(spin_phase_deltas)))
    net_abs_spin_revolutions = abs(net_spin_deg) / 360.0
    total_abs_spin_revolutions = total_abs_spin_deg / 360.0
    signed_spin_rate_hz = (
        net_spin_deg / 360.0 / duration_s
        if np.isfinite(duration_s) and duration_s > 0.0
        else float("nan")
    )
    direction_consistency = (
        net_abs_spin_revolutions / total_abs_spin_revolutions
        if total_abs_spin_revolutions > 0.0
        else 0.0
    )
    helix_to_root_ratio = (
        net_abs_spin_revolutions / net_abs_root_revolutions
        if net_abs_root_revolutions > 0.0
        else float("nan")
    )
    min_spin_fit_r2 = float(np.min(spin_fit_r2))

    if first_fail_category == "none":
        hard_limits = [
            (
                "hook",
                "hook_len_rel_err_max",
                max_hook,
                HELIX_RETENTION_HOOK_REL_ERR_MAX_LIMIT,
            ),
            (
                "hook",
                "local_attach_first_rel_err",
                max_attach,
                HELIX_RETENTION_HOOK_REL_ERR_MAX_LIMIT,
            ),
            (
                "flag",
                "flag_bond_rel_err_max",
                max_bond,
                HELIX_RETENTION_BOND_REL_ERR_MAX_LIMIT,
            ),
            (
                "flag",
                "flag_bend_err_max_deg",
                max_bend,
                HELIX_RETENTION_BEND_ERR_MAX_DEG_LIMIT,
            ),
            (
                "flag",
                "flag_torsion_err_max_deg",
                max_torsion,
                HELIX_RETENTION_TORSION_ERR_MAX_DEG_LIMIT,
            ),
        ]
        for category, name, value, limit in hard_limits:
            if not np.isfinite(value) or value > limit:
                first_fail_step = int(_parse_float(eval_rows[-1].get("step")))
                first_fail_category = category
                first_fail_reason = f"{name}={value:.6g} > {limit:.6g}"
                break

    if first_fail_category == "none" and (
        not np.isfinite(median_abs_spin_rate_hz)
        or median_abs_spin_rate_hz < min_median_abs_spin_rate_hz
    ):
        step_value = _parse_float(eval_rows[-1].get("step"))
        first_fail_step = int(step_value) if np.isfinite(step_value) else len(rows) - 1
        first_fail_category = "motor_no_rotation"
        first_fail_reason = (
            "median_abs_flag_helix_spin_rate_hz="
            f"{median_abs_spin_rate_hz:.6g} < {min_median_abs_spin_rate_hz:.6g}"
        )

    if first_fail_category == "none" and (
        not np.isfinite(net_abs_spin_revolutions)
        or net_abs_spin_revolutions < min_net_abs_spin_revolutions
    ):
        step_value = _parse_float(eval_rows[-1].get("step"))
        first_fail_step = int(step_value) if np.isfinite(step_value) else len(rows) - 1
        first_fail_category = "motor_no_rotation"
        first_fail_reason = (
            "net_abs_flag_helix_spin_revolutions="
            f"{net_abs_spin_revolutions:.6g} < {min_net_abs_spin_revolutions:.6g}"
        )

    if first_fail_category == "none" and (
        not np.isfinite(direction_consistency)
        or direction_consistency < min_direction_consistency
    ):
        step_value = _parse_float(eval_rows[-1].get("step"))
        first_fail_step = int(step_value) if np.isfinite(step_value) else len(rows) - 1
        first_fail_category = "motor_no_rotation"
        first_fail_reason = (
            "flag_helix_spin_direction_consistency="
            f"{direction_consistency:.6g} < {min_direction_consistency:.6g}"
        )

    if first_fail_category == "none" and (
        not np.isfinite(min_spin_fit_r2) or min_spin_fit_r2 < min_helix_spin_fit_r2
    ):
        step_value = _parse_float(eval_rows[-1].get("step"))
        first_fail_step = int(step_value) if np.isfinite(step_value) else len(rows) - 1
        first_fail_category = "motor_spin_metric"
        first_fail_reason = (
            "min_flag_helix_spin_fit_r2="
            f"{min_spin_fit_r2:.6g} < {min_helix_spin_fit_r2:.6g}"
        )

    return {
        "helix_retention_pass": first_fail_category == "none",
        "step_count": len(eval_rows),
        "first_fail_step": first_fail_step,
        "first_fail_category": first_fail_category,
        "first_fail_reason": first_fail_reason,
        "median_abs_flag_phase_rate_hz": median_abs_phase_rate_hz,
        "net_abs_flag_root_revolutions": net_abs_root_revolutions,
        "signed_flag_root_rate_hz": signed_root_rate_hz,
        "median_abs_flag_helix_spin_rate_hz": median_abs_spin_rate_hz,
        "net_abs_flag_helix_spin_revolutions": net_abs_spin_revolutions,
        "signed_flag_helix_spin_rate_hz": signed_spin_rate_hz,
        "flag_helix_spin_direction_consistency": direction_consistency,
        "helix_to_root_net_rotation_ratio": helix_to_root_ratio,
        "min_flag_helix_spin_fit_r2": min_spin_fit_r2,
        "max_hook_len_rel_err": max_hook,
        "max_local_attach_first_rel_err": max_attach,
        "max_flag_bond_rel_err": max_bond,
        "max_flag_bend_err_deg": max_bend,
        "max_flag_torsion_err_deg": max_torsion,
    }
