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
    "final_shape_pass_nonbody",
    "body_shape_pass",
    "shape_pass",
    "first_fail_category_nonbody",
    "body_fail_category",
)
CONTINUOUS_FIELDS = (
    "max_hook_angle_err_deg",
    "max_flag_helix_radius_abs_err_over_b",
    "max_flag_helix_pitch_rel_err",
    "max_net_force_residual_ratio",
    "max_net_torque_residual_ratio",
    "max_motor_force_balance_residual_ratio",
    "max_motor_torque_balance_residual_ratio",
)
MANIFEST_FIELDS = (
    ("kind",),
    ("stage",),
    ("duration_tau",),
    ("dt_star",),
    ("motor_enabled",),
    ("base_config",),
    ("condition_order",),
    ("git", "commit"),
    ("git", "is_clean"),
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
    if math.isnan(left_value) or math.isnan(right_value):
        return math.isnan(left_value) and math.isnan(right_value), None
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
    checks.append(
        _check(
            "summary.identities",
            sorted(left.rows),
            sorted(right.rows),
            set(left.rows) == set(right.rows),
        )
    )
    for identity in sorted(set(left.rows) & set(right.rows)):
        left_row = left.rows[identity]
        right_row = right.rows[identity]
        for field in DISCRETE_FIELDS:
            if field in left_row or field in right_row:
                checks.append(
                    _check(
                        f"{identity}.{field}",
                        left_row.get(field),
                        right_row.get(field),
                        left_row.get(field) == right_row.get(field),
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
