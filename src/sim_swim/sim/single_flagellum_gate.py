"""single flagellum 短時間 motor-on 条件の gate 判定。"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np

REQUIRED_SINGLE_FLAGELLUM_NUMERIC_COLUMNS = (
    "motor_Ta_norm",
    "motor_attach_force_norm",
    "motor_split_residual_norm",
    "flag_root_azimuth_deg",
    "flag_phase_deg",
    "flag_phase_rate_hz",
    "flag_body_phase_diff_deg",
    "local_attach_first_rel_err",
    "local_first_second_rel_err",
    "hook_len_rel_err_max",
    "flag_bond_rel_err_max",
    "flag_bend_err_max_deg",
    "flag_torsion_err_max_deg",
)

REQUIRED_SINGLE_FLAGELLUM_ZERO_COUNTER_COLUMNS = (
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


def _empty_summary(category: str, reason: str) -> dict[str, float | bool | int | str]:
    return {
        "single_flagellum_short_run_pass": False,
        "first_fail_step": -1,
        "first_fail_category": category,
        "first_fail_reason": reason,
        "median_abs_flag_phase_rate_hz": float("nan"),
        "max_hook_len_rel_err": float("nan"),
        "max_local_attach_first_rel_err": float("nan"),
        "max_flag_bond_rel_err": float("nan"),
        "max_flag_bend_err_deg": float("nan"),
        "max_flag_torsion_err_deg": float("nan"),
        "body_shape_pass": False,
        "body_fail_category": "not_evaluated",
    }


def summarize_single_flagellum_short_run(
    rows: Sequence[Mapping[str, str | float | None]],
    body_shape: Mapping[str, str | float | bool | None],
    *,
    skip_initial_steps: int = 1,
    min_median_abs_phase_rate_hz: float = 1.0,
) -> dict[str, float | bool | int | str]:
    """full flagellum の短時間 motor-on 条件を gate 判定する。

    判定は「数値発散しない」「non-body shape gate を通る」「body gate を通る」
    「motor split counter が増えない」「位相変化が実際にある」を同時に要求する。
    """

    if not rows:
        return _empty_summary("empty", "step_summary rows are empty")

    eval_rows = list(rows[int(skip_initial_steps) :])
    if not eval_rows:
        return _empty_summary("empty", "no rows remain after initial-step skip")

    body_shape_pass = _parse_bool(body_shape.get("body_shape_pass"))
    body_fail_category = str(body_shape.get("body_fail_category", "none"))
    if not body_shape_pass:
        summary = _empty_summary("body", body_fail_category)
        summary["body_shape_pass"] = False
        summary["body_fail_category"] = body_fail_category
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
        for col in REQUIRED_SINGLE_FLAGELLUM_NUMERIC_COLUMNS:
            value = _parse_float(row.get(col))
            if not np.isfinite(value):
                first_fail_step = step
                first_fail_category = "nonfinite"
                first_fail_reason = f"non-finite value: {col}"
                break
        if first_fail_category != "none":
            break
        for col in REQUIRED_SINGLE_FLAGELLUM_ZERO_COUNTER_COLUMNS:
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

    phase_rates = np.asarray(
        [_parse_float(row.get("flag_phase_rate_hz")) for row in eval_rows],
        dtype=float,
    )
    median_abs_phase_rate_hz = float(np.median(np.abs(phase_rates)))
    if first_fail_category == "none" and (
        not np.isfinite(median_abs_phase_rate_hz)
        or median_abs_phase_rate_hz < min_median_abs_phase_rate_hz
    ):
        step_value = _parse_float(eval_rows[-1].get("step"))
        first_fail_step = int(step_value) if np.isfinite(step_value) else len(rows) - 1
        first_fail_category = "motor_no_rotation"
        first_fail_reason = (
            "median_abs_flag_phase_rate_hz="
            f"{median_abs_phase_rate_hz:.6g} < {min_median_abs_phase_rate_hz:.6g}"
        )

    def max_col(name: str) -> float:
        values = np.asarray([_parse_float(row.get(name)) for row in eval_rows])
        return float(np.max(values)) if np.isfinite(values).all() else float("nan")

    return {
        "single_flagellum_short_run_pass": first_fail_category == "none",
        "first_fail_step": first_fail_step,
        "first_fail_category": first_fail_category,
        "first_fail_reason": first_fail_reason,
        "median_abs_flag_phase_rate_hz": median_abs_phase_rate_hz,
        "max_hook_len_rel_err": max_col("hook_len_rel_err_max"),
        "max_local_attach_first_rel_err": max_col("local_attach_first_rel_err"),
        "max_flag_bond_rel_err": max_col("flag_bond_rel_err_max"),
        "max_flag_bend_err_deg": max_col("flag_bend_err_max_deg"),
        "max_flag_torsion_err_deg": max_col("flag_torsion_err_max_deg"),
        "body_shape_pass": body_shape_pass,
        "body_fail_category": body_fail_category,
    }
