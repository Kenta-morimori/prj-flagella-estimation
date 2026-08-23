"""Compare bounded Phase 2 qualification artifacts without reading raw CSVs."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[3]


DEFAULT_ATOL = 1.0e-9
DEFAULT_RTOL = 1.0e-6

DISCRETE_FIELDS = (
    "status",
    "expected_steps",
    "completed_steps",
    "completion_pass",
    "finite_pass_all",
    "finite_pass",
    "shape_pass_nonbody",
    "final_shape_pass_nonbody",
    "body_shape_pass",
    "shape_pass",
    "first_fail_category_nonbody",
    "body_fail_category",
    "final_first_fail_category_nonbody",
    "first_fail_category",
)
CONTINUOUS_FIELDS = (
    "max_hook_angle_err_deg",
    "max_flag_helix_radius_abs_err_over_b",
    "max_flag_helix_pitch_rel_err",
    "max_net_force_residual_ratio",
    "max_net_torque_residual_ratio",
    "max_motor_force_balance_residual_ratio",
    "max_motor_torque_balance_residual_ratio",
    "max_hook_len_rel_err",
    "hook_len_rel_err_max",
    "max_flag_bond_rel_err",
    "flag_bond_rel_err_max",
    "flag_bend_err_max_deg",
    "flag_torsion_err_max_deg",
    "body_spring_max_stretch_ratio",
    "body_bend_max_error_deg",
    "body_centerline_max_deviation_um",
    "body_triangle_area_ratio_min",
)
MANIFEST_FIELDS = (
    ("kind",),
    ("stage",),
    ("duration_tau",),
    ("dt_star",),
    ("motor_enabled",),
    ("base_config",),
    ("condition_order",),
)
EXECUTION_ARGUMENT_FIELDS = (
    "mode",
    "duration_s",
    "dt_star",
    "torque_nm",
    "n_flagella",
    "attach_seed",
    "phase_seed",
)
SUCCESS_TRUE_FIELDS = frozenset(
    {
        "completion_pass",
        "finite_pass_all",
        "finite_pass",
        "shape_pass_nonbody",
        "final_shape_pass_nonbody",
        "body_shape_pass",
        "shape_pass",
    }
)
FAILURE_CATEGORY_FIELDS = frozenset(
    {
        "first_fail_category_nonbody",
        "body_fail_category",
        "final_first_fail_category_nonbody",
        "first_fail_category",
    }
)


@dataclass(frozen=True)
class Campaign:
    root: Path
    manifest_path: Path
    manifest: dict[str, Any]
    summary_path: Path
    rows: dict[str, dict[str, str]]


def _read_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"JSON object required: {path}")
    return value


def _manifest_under(path: Path) -> Path:
    if path.is_file():
        if path.name != "run_manifest.json":
            raise ValueError(f"Expected run_manifest.json, got: {path}")
        return path
    direct = path / "run_manifest.json"
    if direct.is_file():
        return direct
    candidates = sorted(path.glob("**/run_manifest.json"))
    if len(candidates) != 1:
        raise ValueError(
            f"Expected exactly one run_manifest.json under {path}; found {len(candidates)}"
        )
    return candidates[0]


def load_campaign(path: Path) -> Campaign:
    """Load one campaign through its compact summary and manifest."""

    manifest_path = _manifest_under(path.resolve())
    root = manifest_path.parent
    manifest = _read_json(manifest_path)
    raw_summary = manifest.get("summary_csv")
    candidates = []
    if isinstance(raw_summary, str):
        candidates.append(Path(raw_summary))
    candidates.append(root / "summary.csv")
    summary_path = next(
        (candidate for candidate in candidates if candidate.is_file()), None
    )
    if summary_path is None:
        raise ValueError(f"summary.csv is missing for campaign: {root}")
    with summary_path.open(encoding="utf-8", newline="") as handle:
        raw_rows = list(csv.DictReader(handle))
    rows: dict[str, dict[str, str]] = {}
    for row in raw_rows:
        key = row.get("profile") or row.get("condition_id")
        if not key:
            raise ValueError(f"summary row lacks profile/condition_id: {summary_path}")
        if key in rows:
            raise ValueError(f"duplicate summary identity {key!r}: {summary_path}")
        rows[key] = row
    return Campaign(root, manifest_path, manifest, summary_path, rows)


def _nested(data: dict[str, Any], path: tuple[str, ...]) -> Any:
    value: Any = data
    for key in path:
        if not isinstance(value, dict):
            return None
        value = value.get(key)
    return value


def _as_float(value: Any) -> float | None:
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed


def _numeric_match(
    left: Any, right: Any, *, atol: float, rtol: float
) -> tuple[bool, float | None]:
    left_value = _as_float(left)
    right_value = _as_float(right)
    if left_value is None or right_value is None:
        return left == right, None
    if not math.isfinite(left_value) or not math.isfinite(right_value):
        return False, None
    difference = abs(left_value - right_value)
    allowed = max(atol, rtol * max(abs(left_value), abs(right_value)))
    return difference <= allowed, difference


def _check(
    name: str, left: Any, right: Any, passed: bool, **extra: Any
) -> dict[str, Any]:
    return {
        "name": name,
        "left": left,
        "right": right,
        "status": "pass" if passed else "fail",
        **extra,
    }


def _is_true(value: Any) -> bool:
    return value is True or value == "True"


def _is_no_failure(value: Any) -> bool:
    return value in (None, "", "none")


def _discrete_pass(field: str, left: Any, right: Any) -> bool:
    """Require both equality and a success state for qualification fields."""

    if left != right:
        return False
    if field == "status":
        return left == "completed"
    if field in SUCCESS_TRUE_FIELDS:
        return left in (None, "") or _is_true(left)
    if field in FAILURE_CATEGORY_FIELDS:
        return _is_no_failure(left)
    return True


def _config_sha256(config: str) -> str | None:
    path = (REPO_ROOT / config).resolve()
    try:
        path.relative_to((REPO_ROOT / "conf" / "phase2_sweeps").resolve())
    except ValueError:
        return None
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def _parallel_output_path(
    record: dict[str, Any], manifest: dict[str, Any], job_manifest_path: Path
) -> Path | None:
    raw = record.get("output_dir")
    if not isinstance(raw, str):
        return None
    direct = Path(raw)
    if direct.is_dir():
        return direct
    root = manifest.get("output_root")
    if not isinstance(root, str):
        return direct
    try:
        relative = direct.relative_to(Path(root))
    except ValueError:
        return direct
    return job_manifest_path.parent / relative


def _run_summary_for(
    campaign: Campaign, identity: str, row: dict[str, str]
) -> Path | None:
    candidates = []
    raw_output_dir = row.get("output_dir")
    if raw_output_dir:
        candidates.append(Path(raw_output_dir) / "run_summary.json")
    for condition in campaign.manifest.get("conditions", []):
        if not isinstance(condition, dict) or condition.get("condition_id") != identity:
            continue
        raw_path = condition.get("run_summary_json")
        if isinstance(raw_path, str):
            candidates.append(Path(raw_path))
    candidates.append(campaign.root / identity / "run_summary.json")
    return next((candidate for candidate in candidates if candidate.is_file()), None)


def _execution_checks(campaign: Campaign, identity: str) -> list[dict[str, Any]]:
    path = _run_summary_for(campaign, identity, campaign.rows[identity])
    prefix = f"{identity}.run_summary"
    if path is None:
        return [_check(prefix, "missing", "present", False)]
    execution = _read_json(path).get("execution")
    if not isinstance(execution, dict):
        return [_check(prefix + ".execution", "missing", "present", False)]
    expected_steps = execution.get("expected_total_steps")
    row_count = execution.get("row_count")
    return [
        _check(
            prefix + ".status",
            execution.get("status"),
            "completed",
            execution.get("status") == "completed",
        ),
        _check(
            prefix + ".step_indices_contiguous_from_zero",
            execution.get("step_indices_contiguous_from_zero"),
            True,
            execution.get("step_indices_contiguous_from_zero") is True,
        ),
        _check(
            prefix + ".expected_total_steps",
            expected_steps,
            row_count,
            isinstance(expected_steps, int) and expected_steps == row_count,
        ),
    ]


def _condition_order(manifest: dict[str, Any]) -> list[str] | None:
    """Return a valid, non-empty manifest condition order when present."""

    raw = manifest.get("condition_order")
    if (
        not isinstance(raw, list)
        or not raw
        or any(not isinstance(identity, str) or not identity for identity in raw)
        or len(set(raw)) != len(raw)
    ):
        return None
    return raw


def compare_campaigns(
    left_path: Path,
    right_path: Path,
    *,
    atol: float = DEFAULT_ATOL,
    rtol: float = DEFAULT_RTOL,
) -> dict[str, Any]:
    """Compare two serial-like campaign outputs using their compact artifacts."""

    if atol < 0 or rtol < 0:
        raise ValueError("atol and rtol must be non-negative")
    left = load_campaign(left_path)
    right = load_campaign(right_path)
    checks: list[dict[str, Any]] = []
    for field in MANIFEST_FIELDS:
        left_value = _nested(left.manifest, field)
        right_value = _nested(right.manifest, field)
        checks.append(
            _check(
                "manifest." + ".".join(field),
                left_value,
                right_value,
                left_value == right_value,
            )
        )
    left_commit = _nested(left.manifest, ("git", "commit"))
    right_commit = _nested(right.manifest, ("git", "commit"))
    checks.append(
        _check(
            "manifest.git.commit",
            left_commit,
            right_commit,
            isinstance(left_commit, str)
            and bool(left_commit)
            and left_commit == right_commit,
        )
    )
    left_clean = _nested(left.manifest, ("git", "is_clean"))
    right_clean = _nested(right.manifest, ("git", "is_clean"))
    checks.append(
        _check(
            "manifest.git.is_clean",
            left_clean,
            right_clean,
            left_clean is True and right_clean is True,
        )
    )
    for field in EXECUTION_ARGUMENT_FIELDS:
        left_value = _nested(left.manifest, ("args", field))
        right_value = _nested(right.manifest, ("args", field))
        if left_value is not None or right_value is not None:
            checks.append(
                _check(
                    "manifest.args." + field,
                    left_value,
                    right_value,
                    left_value == right_value,
                )
            )
    checks.append(
        _check(
            "summary.identities",
            sorted(left.rows),
            sorted(right.rows),
            set(left.rows) == set(right.rows),
        )
    )
    for side, campaign in (("left", left), ("right", right)):
        expected = _condition_order(campaign.manifest)
        actual = sorted(campaign.rows)
        checks.append(
            _check(
                f"{side}.manifest.condition_order",
                campaign.manifest.get("condition_order"),
                "non-empty unique condition identities",
                expected is not None,
            )
        )
        checks.append(
            _check(
                f"{side}.summary.identities",
                actual,
                sorted(expected) if expected is not None else [],
                expected is not None and actual == sorted(expected),
            )
        )
    for identity in sorted(set(left.rows) & set(right.rows)):
        if left.manifest.get("kind") == "shape_stability_grid":
            checks.extend(
                {**check, "name": "left." + check["name"]}
                for check in _execution_checks(left, identity)
            )
        if right.manifest.get("kind") == "shape_stability_grid":
            checks.extend(
                {**check, "name": "right." + check["name"]}
                for check in _execution_checks(right, identity)
            )
        left_row = left.rows[identity]
        right_row = right.rows[identity]
        for field in DISCRETE_FIELDS:
            if field in left_row or field in right_row:
                checks.append(
                    _check(
                        f"{identity}.{field}",
                        left_row.get(field),
                        right_row.get(field),
                        _discrete_pass(
                            field, left_row.get(field), right_row.get(field)
                        ),
                    )
                )
        for field in CONTINUOUS_FIELDS:
            if field in left_row or field in right_row:
                passed, difference = _numeric_match(
                    left_row.get(field), right_row.get(field), atol=atol, rtol=rtol
                )
                checks.append(
                    _check(
                        f"{identity}.{field}",
                        left_row.get(field),
                        right_row.get(field),
                        passed,
                        absolute_difference=difference,
                        atol=atol,
                        rtol=rtol,
                    )
                )
    return {
        "kind": "phase2_qualification_comparison",
        "left": {
            "root": str(left.root),
            "manifest": str(left.manifest_path),
            "summary": str(left.summary_path),
        },
        "right": {
            "root": str(right.root),
            "manifest": str(right.manifest_path),
            "summary": str(right.summary_path),
        },
        "atol": atol,
        "rtol": rtol,
        "checks": checks,
        "status": "PASS"
        if all(item["status"] == "pass" for item in checks)
        else "FAIL",
    }


def compare_parallel_job(
    job_manifest_path: Path,
    serial_campaigns: dict[str, Path],
    *,
    atol: float = DEFAULT_ATOL,
    rtol: float = DEFAULT_RTOL,
) -> dict[str, Any]:
    """Compare every completed child in a parallel job with its serial reference."""

    manifest = _read_json(job_manifest_path.resolve())
    checks: list[dict[str, Any]] = []
    checks.append(
        _check(
            "job.status",
            manifest.get("status"),
            "succeeded",
            manifest.get("status") == "succeeded",
        )
    )
    checks.append(
        _check(
            "job.failed_configs",
            manifest.get("failed_configs"),
            [],
            manifest.get("failed_configs") == [],
        )
    )
    records = manifest.get("configs")
    if not isinstance(records, list):
        raise ValueError(f"job manifest configs must be a list: {job_manifest_path}")
    seen_configs: set[str] = set()
    for record in records:
        if not isinstance(record, dict):
            raise ValueError(
                f"job manifest config record must be an object: {job_manifest_path}"
            )
        config = record.get("config")
        if not isinstance(config, str):
            raise ValueError(
                f"job manifest config path is missing: {job_manifest_path}"
            )
        seen_configs.add(config)
        config_sha = record.get("config_sha256")
        if config_sha is not None:
            checks.append(
                _check(
                    f"{config}.config_sha256",
                    config_sha,
                    _config_sha256(config),
                    config_sha == _config_sha256(config),
                )
            )
        checks.append(
            _check(
                f"{config}.launcher_status",
                record.get("status"),
                "succeeded",
                record.get("status") == "succeeded",
            )
        )
        reference = serial_campaigns.get(config)
        if reference is None:
            checks.append(_check(f"{config}.serial_reference", None, "provided", False))
            continue
        output_dir = _parallel_output_path(record, manifest, job_manifest_path)
        if output_dir is None:
            checks.append(
                _check(f"{config}.parallel_output_dir", None, "provided", False)
            )
            continue
        pair = compare_campaigns(reference, output_dir, atol=atol, rtol=rtol)
        for check in pair["checks"]:
            checks.append({**check, "name": f"{config}.{check['name']}"})
    for config in sorted(set(serial_campaigns) - seen_configs):
        checks.append(_check(f"{config}.parallel_record", "missing", "present", False))
    return {
        "kind": "phase2_parallel_qualification_comparison",
        "job_manifest": str(job_manifest_path.resolve()),
        "serial_campaigns": {
            key: str(value) for key, value in serial_campaigns.items()
        },
        "atol": atol,
        "rtol": rtol,
        "checks": checks,
        "status": "PASS"
        if all(item["status"] == "pass" for item in checks)
        else "FAIL",
    }


def write_report(report: dict[str, Any], output_dir: Path) -> tuple[Path, Path]:
    """Write machine-readable JSON and a compact per-check CSV report."""

    output_dir.mkdir(parents=True, exist_ok=False)
    json_path = output_dir / "qualification_report.json"
    csv_path = output_dir / "qualification_report.csv"
    json_path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    with csv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["name", "status", "left", "right", "absolute_difference"],
        )
        writer.writeheader()
        for check in report["checks"]:
            writer.writerow({key: check.get(key, "") for key in writer.fieldnames})
    (output_dir / "manifest.json").write_text(
        json.dumps(
            {
                "kind": "phase2_qualification_report",
                "status": report["status"],
                "report_json": str(json_path),
                "report_csv": str(csv_path),
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (output_dir / "run.log").write_text(
        f"phase2 qualification comparison {report['status']}\n", encoding="utf-8"
    )
    return json_path, csv_path
