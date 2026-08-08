"""Phase 2 run diagnostics を軽量な run_summary.json へ集約する。"""

from __future__ import annotations

import csv
import json
import math
from collections import Counter
from datetime import datetime
from pathlib import Path
from typing import Any, Mapping, Sequence
from zoneinfo import ZoneInfo

from sim_swim.sim.body_shape_gate import (
    BODY_BEND_ERR_MAX_DEG_LIMIT,
    BODY_CENTERLINE_DEV_UM_MAX_LIMIT,
    BODY_SPRING_STRETCH_RATIO_MAX_LIMIT,
    BODY_TRIANGLE_AREA_RATIO_MIN_LIMIT,
)

SCHEMA_VERSION = "1.0"
EPISODE_STORAGE_LIMIT = 32
PERSISTENT_MIN_CONSECUTIVE_FAIL_SAMPLES = 3
MAX_RUN_SUMMARY_BYTES = 64 * 1024


def _float(value: object) -> float:
    try:
        return float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return float("nan")


def _bool(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def resolve_run_dir(input_dir: Path) -> Path:
    """``step_summary.csv`` を含む run directory を解決する。"""
    input_dir = Path(input_dir)
    if (input_dir / "step_summary.csv").is_file():
        return input_dir
    if (input_dir / "sim" / "step_summary.csv").is_file():
        return input_dir / "sim"
    raise FileNotFoundError(f"step_summary.csv is missing below: {input_dir}")


def _read_json(path: Path) -> dict[str, Any] | None:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None
    return data if isinstance(data, dict) else None


def _manifest_time(run_dir: Path) -> dict[str, Any] | None:
    for candidate in (run_dir / "manifest.json", run_dir.parent / "manifest.json"):
        manifest = _read_json(candidate)
        if isinstance(manifest, dict) and isinstance(manifest.get("time"), dict):
            return dict(manifest["time"])
    campaign = _read_json(run_dir.parent / "run_manifest.json")
    if not isinstance(campaign, dict):
        return None
    for condition in campaign.get("conditions", []):
        if (
            isinstance(condition, dict)
            and Path(str(condition.get("output_dir", ""))).name == run_dir.name
        ):
            if isinstance(condition.get("time"), dict):
                return dict(condition["time"])
    return dict(campaign["time"]) if isinstance(campaign.get("time"), dict) else None


def _execution(
    rows: Sequence[Mapping[str, str]], time: Mapping[str, Any] | None
) -> dict[str, Any]:
    steps = [_float(row.get("step")) for row in rows]
    times = [_float(row.get("t_s")) for row in rows]
    contiguous = bool(rows) and all(
        math.isfinite(step) and int(step) == index for index, step in enumerate(steps)
    )
    result: dict[str, Any] = {
        "status": "unknown",
        "row_count": len(rows),
        "step_indices_contiguous_from_zero": contiguous,
        "observed_first_t_s": times[0] if times and math.isfinite(times[0]) else None,
        "observed_final_t_s": times[-1] if times and math.isfinite(times[-1]) else None,
    }
    if not time:
        result["reason"] = "time manifest is unavailable"
        return result
    expected_steps = _float(time.get("total_steps"))
    expected_final = _float(time.get("final_step_summary_t_s"))
    result["expected_total_steps"] = (
        int(expected_steps) if math.isfinite(expected_steps) else None
    )
    result["expected_final_step_summary_t_s"] = (
        expected_final if math.isfinite(expected_final) else None
    )
    if not math.isfinite(expected_steps) or not math.isfinite(expected_final):
        result["reason"] = "time manifest lacks total_steps or final_step_summary_t_s"
        return result
    completed = (
        contiguous
        and len(rows) == int(expected_steps)
        and bool(times)
        and math.isfinite(times[-1])
        and math.isclose(times[-1], expected_final, rel_tol=1e-9, abs_tol=1e-12)
    )
    result["status"] = "completed" if completed else "partial"
    result["reason"] = (
        "matched time manifest"
        if completed
        else "observed rows do not match time manifest"
    )
    return result


def _sample_spacing(times: Sequence[float]) -> dict[str, float | None]:
    gaps = sorted(
        later - earlier for earlier, later in zip(times, times[1:]) if later > earlier
    )
    if not gaps:
        return {"min_s": None, "median_s": None, "max_s": None}
    midpoint = len(gaps) // 2
    median = (
        gaps[midpoint] if len(gaps) % 2 else (gaps[midpoint - 1] + gaps[midpoint]) / 2
    )
    return {"min_s": gaps[0], "median_s": median, "max_s": gaps[-1]}


def _episodes(
    rows: Sequence[Mapping[str, str]],
    passed: Sequence[bool],
    *,
    category_column: str | None = None,
) -> tuple[list[dict[str, Any]], int]:
    times = [_float(row.get("t_s")) for row in rows]
    spans: list[tuple[int, int]] = []
    start: int | None = None
    for index, is_pass in enumerate(passed):
        if not is_pass and start is None:
            start = index
        elif is_pass and start is not None:
            spans.append((start, index))
            start = None
    if start is not None:
        spans.append((start, len(rows)))
    episodes: list[dict[str, Any]] = []
    for begin, end in spans:
        begin_t = times[begin]
        end_t = times[end - 1]
        categories = Counter(
            str(row.get(category_column, "unspecified"))[:128]
            for row in rows[begin:end]
            if category_column and str(row.get(category_column, ""))
        )
        episodes.append(
            {
                "start_step": int(_float(rows[begin].get("step"))),
                "last_observed_fail_step": int(_float(rows[end - 1].get("step"))),
                "start_t_s": begin_t if math.isfinite(begin_t) else None,
                "end_t_s": end_t if math.isfinite(end_t) else None,
                "next_observed_pass_t_s": (
                    times[end]
                    if end < len(times) and math.isfinite(times[end])
                    else None
                ),
                "observed_fail_sample_count": end - begin,
                "observed_duration_s": (
                    max(0.0, end_t - begin_t)
                    if math.isfinite(begin_t) and math.isfinite(end_t)
                    else None
                ),
                "persistent_observed": (end - begin)
                >= PERSISTENT_MIN_CONSECUTIVE_FAIL_SAMPLES,
                "category_counts": dict(sorted(categories.items())),
            }
        )
    if len(episodes) <= EPISODE_STORAGE_LIMIT:
        return episodes, 0
    ranked = sorted(
        range(len(episodes)),
        key=lambda index: (
            -float(episodes[index]["observed_duration_s"] or 0.0),
            index,
        ),
    )
    kept = set(range(8)) | set(range(len(episodes) - 8, len(episodes)))
    for index in ranked:
        if len(kept) >= EPISODE_STORAGE_LIMIT:
            break
        kept.add(index)
    selected = [episode for index, episode in enumerate(episodes) if index in kept]
    return selected, len(episodes) - len(selected)


def _gate_summary(
    rows: Sequence[Mapping[str, str]],
    passed: Sequence[bool] | None,
    *,
    category_column: str | None = None,
    unavailable_reason: str | None = None,
) -> dict[str, Any]:
    if passed is None:
        return {"status": "unavailable", "reason": unavailable_reason}
    episodes, omitted = _episodes(rows, passed, category_column=category_column)
    failed = [index for index, value in enumerate(passed) if not value]
    return {
        "status": "available",
        "final_pass": bool(passed[-1]) if passed else None,
        "any_fail": bool(failed),
        "first_observed_fail_t_s": (
            _float(rows[failed[0]].get("t_s")) if failed else None
        ),
        "last_observed_fail_t_s": (
            _float(rows[failed[-1]].get("t_s")) if failed else None
        ),
        "observed_fail_sample_count": len(failed),
        "episode_count": len(episodes) + omitted,
        "stored_episode_count": len(episodes),
        "omitted_episode_count": omitted,
        "recovered_after_fail": bool(failed) and bool(passed[-1]),
        "episodes": episodes,
    }


def _step_gate(
    rows: Sequence[Mapping[str, str]],
    column: str,
    *,
    category_column: str | None = None,
) -> dict[str, Any]:
    if not rows:
        return _gate_summary(rows, None, unavailable_reason="step_summary.csv is empty")
    if any(column not in row for row in rows):
        return _gate_summary(
            rows,
            None,
            unavailable_reason=f"step_summary.csv lacks required column: {column}",
        )
    return _gate_summary(
        rows,
        [_bool(row.get(column)) for row in rows],
        category_column=category_column,
    )


def _body_gate(rows: Sequence[Mapping[str, str]]) -> dict[str, Any]:
    if not rows:
        return _gate_summary(
            rows,
            None,
            unavailable_reason="body_constraint_diagnostics.csv is missing or empty",
        )
    initial_area = _float(rows[0].get("body_triangle_area_min"))
    passed: list[bool] = []
    normalized: list[dict[str, str]] = []
    for row in rows:
        stretch = _float(row.get("body_spring_max_stretch_ratio"))
        bend = _float(row.get("body_bend_max_error_deg"))
        centerline = _float(row.get("body_centerline_max_deviation_um"))
        area = _float(row.get("body_triangle_area_min"))
        if (
            not all(
                math.isfinite(value)
                for value in (stretch, bend, centerline, area, initial_area)
            )
            or initial_area <= 0.0
        ):
            category = "body_nonfinite"
        elif stretch > BODY_SPRING_STRETCH_RATIO_MAX_LIMIT:
            category = "body_spring"
        elif bend > BODY_BEND_ERR_MAX_DEG_LIMIT:
            category = "body_bend"
        elif centerline > BODY_CENTERLINE_DEV_UM_MAX_LIMIT:
            category = "body_centerline"
        elif area / max(initial_area, 1e-30) < BODY_TRIANGLE_AREA_RATIO_MIN_LIMIT:
            category = "body_area"
        else:
            category = "none"
        copy = dict(row)
        copy["body_fail_category"] = category
        normalized.append(copy)
        passed.append(category == "none")
    return _gate_summary(normalized, passed, category_column="body_fail_category")


def _extrema(rows: Sequence[Mapping[str, str]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for metric in (
        "hook_len_rel_err_max",
        "flag_bond_rel_err_max",
        "flag_bend_err_max_deg",
        "flag_torsion_err_max_deg",
    ):
        candidates = [
            (index, _float(row.get(metric))) for index, row in enumerate(rows)
        ]
        finite = [(index, value) for index, value in candidates if math.isfinite(value)]
        if not finite:
            result[metric] = {"value": None, "t_s": None}
            continue
        index, value = max(finite, key=lambda item: item[1])
        result[metric] = {"value": value, "t_s": _float(rows[index].get("t_s"))}
    return result


def build_run_summary(
    input_dir: Path, *, time_manifest: Mapping[str, Any] | None = None
) -> dict[str, Any]:
    """既存 diagnostics から、時系列を複製しない run summary を構築する。"""
    run_dir = resolve_run_dir(input_dir)
    step_path = run_dir / "step_summary.csv"
    step_rows = _rows(step_path)
    body_path = run_dir / "body_constraint_diagnostics.csv"
    body_rows = _rows(body_path) if body_path.is_file() else []
    effective_time = (
        dict(time_manifest) if time_manifest is not None else _manifest_time(run_dir)
    )
    return {
        "schema_version": SCHEMA_VERSION,
        "kind": "phase2_run_summary",
        "created_at": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(),
        "input": {
            "run_dir": str(run_dir),
            "step_summary_csv": str(step_path),
            "body_constraint_diagnostics_csv": (
                str(body_path) if body_path.is_file() else None
            ),
        },
        "execution": _execution(step_rows, effective_time),
        "sampling": {
            "step_summary_row_count": len(step_rows),
            "time_spacing_s": _sample_spacing(
                [_float(row.get("t_s")) for row in step_rows]
            ),
            "episode_definition": "maximal consecutive observed fail samples",
            "persistent_observed_min_consecutive_fail_samples": PERSISTENT_MIN_CONSECUTIVE_FAIL_SAMPLES,
            "episode_storage_limit_per_gate": EPISODE_STORAGE_LIMIT,
        },
        "gates": {
            "finite": _step_gate(step_rows, "finite_pass"),
            "shape_nonbody": _step_gate(
                step_rows,
                "shape_pass_nonbody",
                category_column="first_fail_category_nonbody",
            ),
            "shape_body": _body_gate(body_rows),
        },
        "extrema": _extrema(step_rows),
        "source_file_size_bytes": {
            "step_summary_csv": step_path.stat().st_size,
            "body_constraint_diagnostics_csv": (
                body_path.stat().st_size if body_path.is_file() else None
            ),
        },
    }


def write_run_summary(
    input_dir: Path,
    *,
    time_manifest: Mapping[str, Any] | None = None,
    overwrite: bool = False,
) -> Path:
    run_dir = resolve_run_dir(input_dir)
    output = run_dir / "run_summary.json"
    if output.exists() and not overwrite:
        raise FileExistsError(
            f"run_summary.json already exists: {output}; use --overwrite"
        )
    text = (
        json.dumps(
            build_run_summary(run_dir, time_manifest=time_manifest),
            ensure_ascii=False,
            indent=2,
        )
        + "\n"
    )
    encoded = text.encode("utf-8")
    if len(encoded) > MAX_RUN_SUMMARY_BYTES:
        raise ValueError(
            f"run_summary.json would exceed {MAX_RUN_SUMMARY_BYTES} bytes: {output}"
        )
    output.write_bytes(encoded)
    return output
