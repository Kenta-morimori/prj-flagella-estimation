"""body-only / body-shape 診断の pass/fail 判定。"""

from __future__ import annotations

import csv
from collections.abc import Mapping, Sequence
from pathlib import Path

import numpy as np

BODY_SPRING_STRETCH_RATIO_MAX_LIMIT = 1.0
BODY_BEND_ERR_MAX_DEG_LIMIT = 60.0
BODY_CENTERLINE_DEV_UM_MAX_LIMIT = 2.0
BODY_TRIANGLE_AREA_RATIO_MIN_LIMIT = 0.5


def _parse_float(value: str | float | None) -> float:
    if value is None:
        return float("nan")
    if isinstance(value, float):
        return value
    text = str(value).strip()
    if not text:
        return float("nan")
    try:
        return float(text)
    except ValueError:
        return float("nan")


def _empty_result(category: str) -> dict[str, float | bool | str]:
    return {
        "body_shape_pass": False,
        "body_fail_category": category,
        "body_spring_max_stretch_ratio": float("nan"),
        "body_bend_max_error_deg": float("nan"),
        "body_centerline_max_deviation_um": float("nan"),
        "body_triangle_area_ratio_min": float("nan"),
    }


def summarize_body_shape_diagnostics(
    rows: Sequence[Mapping[str, str | float | None]],
) -> dict[str, float | bool | str]:
    """body diagnostics rows から body shape gate の結果を返す。

    Args:
        rows: `body_constraint_diagnostics.csv` を読み込んだ行列。

    Returns:
        `body_shape_pass`, `body_fail_category`, 主要最大値を含む辞書。
    """

    if not rows:
        return _empty_result("body_diag_empty")

    stretch_vals = [_parse_float(r.get("body_spring_max_stretch_ratio")) for r in rows]
    bend_vals = [_parse_float(r.get("body_bend_max_error_deg")) for r in rows]
    centerline_vals = [
        _parse_float(r.get("body_centerline_max_deviation_um")) for r in rows
    ]
    area_min_vals = [_parse_float(r.get("body_triangle_area_min")) for r in rows]

    spring_max = max(stretch_vals)
    bend_max = max(bend_vals)
    centerline_max = max(centerline_vals)
    init_area_min = area_min_vals[0]
    min_area_min = min(area_min_vals)
    area_ratio_min = (
        min_area_min / max(init_area_min, 1e-30)
        if init_area_min > 0.0
        else float("nan")
    )

    if not all(
        np.isfinite(v) for v in [spring_max, bend_max, centerline_max, area_ratio_min]
    ):
        return {
            "body_shape_pass": False,
            "body_fail_category": "body_nonfinite",
            "body_spring_max_stretch_ratio": spring_max,
            "body_bend_max_error_deg": bend_max,
            "body_centerline_max_deviation_um": centerline_max,
            "body_triangle_area_ratio_min": area_ratio_min,
        }

    if spring_max > BODY_SPRING_STRETCH_RATIO_MAX_LIMIT:
        category = "body_spring"
        passed = False
    elif bend_max > BODY_BEND_ERR_MAX_DEG_LIMIT:
        category = "body_bend"
        passed = False
    elif centerline_max > BODY_CENTERLINE_DEV_UM_MAX_LIMIT:
        category = "body_centerline"
        passed = False
    elif area_ratio_min < BODY_TRIANGLE_AREA_RATIO_MIN_LIMIT:
        category = "body_area"
        passed = False
    else:
        category = "none"
        passed = True

    return {
        "body_shape_pass": passed,
        "body_fail_category": category,
        "body_spring_max_stretch_ratio": spring_max,
        "body_bend_max_error_deg": bend_max,
        "body_centerline_max_deviation_um": centerline_max,
        "body_triangle_area_ratio_min": area_ratio_min,
    }


def summarize_body_shape_diagnostics_csv(
    csv_path: Path,
) -> dict[str, float | bool | str]:
    """body diagnostics CSV から body shape gate の結果を返す。

    Args:
        csv_path: `body_constraint_diagnostics.csv` のパス。

    Returns:
        `summarize_body_shape_diagnostics` と同じ辞書。
    """

    if not csv_path.is_file():
        return _empty_result("body_diag_missing")
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        return summarize_body_shape_diagnostics(list(csv.DictReader(handle)))
