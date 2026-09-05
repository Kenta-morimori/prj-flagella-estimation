"""Assess the Issue #61 2015-project 1τ tracking-reference campaign."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


TORQUES_NM = (1.0e-21, 2.5e-20, 1.0e-19)
EXPECTED_DT_STAR = 1.0e-5
EXPECTED_DURATION_TAU = 1.0


def _read_json(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"JSON object required: {path}")
    return data


def _float(value: Any) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float("nan")


def _close(actual: Any, expected: float) -> bool:
    return math.isclose(_float(actual), expected, rel_tol=1.0e-12, abs_tol=1.0e-30)


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _validate_manifest(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    if manifest.get("kind") != "stage_a_2015" or manifest.get("issue") != 61:
        raise ValueError("run root is not an Issue #61 2015 Stage A campaign")
    for name, expected in (
        ("stage", "motor_on"),
        ("duration_tau", EXPECTED_DURATION_TAU),
        ("dt_star", EXPECTED_DT_STAR),
    ):
        actual = manifest.get(name)
        if isinstance(expected, str):
            if actual != expected:
                raise ValueError(
                    f"run manifest {name}={actual!r}; expected {expected!r}"
                )
        elif not _close(actual, expected):
            raise ValueError(f"run manifest {name}={actual!r}; expected {expected}")
    if manifest.get("profiles") not in (None, ["project"]):
        raise ValueError("Issue #61 requires the 2015 project profile only")
    if manifest.get("link_reference_torque") is not True:
        raise ValueError("Issue #61 requires link_reference_torque=true")
    evidence = manifest.get("reference_evidence")
    if not isinstance(evidence, list):
        raise ValueError("run manifest lacks reference_evidence")
    for item in evidence:
        if not isinstance(item, dict) or not isinstance(
            item.get("manifest_sha256"), str
        ):
            raise ValueError("reference_evidence entry lacks manifest_sha256")
    conditions = manifest.get("conditions")
    if not isinstance(conditions, list) or len(conditions) != len(TORQUES_NM):
        raise ValueError("Issue #61 requires exactly three condition manifests")
    torques = {
        _float(item.get("motor_torque_Nm"))
        for item in conditions
        if isinstance(item, dict)
    }
    if torques != set(TORQUES_NM):
        raise ValueError(f"unexpected torque grid: {sorted(torques)}")
    for condition in conditions:
        if not isinstance(condition, dict):
            raise ValueError("condition manifest must be an object")
        overrides = condition.get("config_overrides")
        if not isinstance(overrides, dict):
            raise ValueError("condition manifest lacks config_overrides")
        torque = condition.get("motor_torque_Nm")
        if not _close(overrides.get("motor.torque_Nm"), _float(torque)):
            raise ValueError("motor torque override disagrees with condition")
        if not _close(overrides.get("motor.reference_torque_Nm"), _float(torque)):
            raise ValueError("reference torque must equal motor torque")
        if overrides.get("time.scale_policy") != "reference_torque":
            raise ValueError("condition lacks reference_torque time policy")
    return [dict(item) for item in conditions]


def _first_gate_failure(summary: dict[str, Any]) -> dict[str, Any] | None:
    execution = summary.get("execution", {})
    if execution.get("status") != "completed":
        return {"criterion": "execution", "t_s": None, "step": None}
    gates = summary.get("gates", {})
    for name in ("finite", "shape_nonbody", "shape_body"):
        gate = gates.get(name, {})
        if gate.get("status") != "available" or gate.get("any_fail"):
            return {
                "criterion": name,
                "t_s": gate.get("first_observed_fail_t_s"),
                "step": None,
            }
    return None


def _threshold_failures(row: dict[str, str], thresholds: dict[str, Any]) -> list[str]:
    failures: list[str] = []
    if (
        row.get("status") != "completed"
        or row.get("completion_pass", "").lower() != "true"
    ):
        failures.append("completion")
    if row.get("finite_pass_all", "").lower() != "true":
        failures.append("finite")
    for metric, raw_limit in thresholds.items():
        value, limit = _float(row.get(metric)), _float(raw_limit)
        if not math.isfinite(value) or value > limit:
            failures.append(metric)
    return failures


def analyze(*, run_root: Path, threshold_contract: Path, output_dir: Path) -> Path:
    """Write a compact PASS/FAIL decision without restarting simulations."""
    manifest = _read_json(run_root / "run_manifest.json")
    conditions = _validate_manifest(manifest)
    threshold_data = _read_json_or_yaml(threshold_contract)
    thresholds = threshold_data.get("thresholds")
    if threshold_data.get("status") != "locked" or not isinstance(thresholds, dict):
        raise ValueError("threshold contract must be locked and contain thresholds")
    summary_by_id = {
        row.get("condition_id"): row for row in _rows(run_root / "summary.csv")
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, Any]] = []
    for condition in sorted(
        conditions, key=lambda item: _float(item["motor_torque_Nm"])
    ):
        condition_id = str(condition["condition_id"])
        row = summary_by_id.get(condition_id)
        if row is None:
            raise ValueError(f"summary.csv is missing {condition_id}")
        condition_dir = Path(str(condition["output_dir"]))
        summary_path = condition_dir / "run_summary.json"
        gate_failure = (
            _first_gate_failure(_read_json(summary_path))
            if summary_path.is_file()
            else {"criterion": "run_summary_missing", "t_s": None, "step": None}
        )
        threshold_failures = _threshold_failures(row, thresholds)
        failures = (
            [gate_failure["criterion"]] if gate_failure else []
        ) + threshold_failures
        first = gate_failure or (
            {"criterion": threshold_failures[0], "t_s": None, "step": None}
            if threshold_failures
            else None
        )
        records.append(
            {
                "condition_id": condition_id,
                "motor_torque_Nm": _float(condition["motor_torque_Nm"]),
                "tau_s": condition.get("time", {}).get("tau_s"),
                "dt_internal_s": condition.get("time", {}).get("dt_internal_s"),
                "total_steps": condition.get("time", {}).get("total_steps"),
                "wall_time_s": _float(row.get("wall_time_s")),
                "steps_per_s": _float(row.get("steps_per_s")),
                "strict_pass": not failures,
                "first_failing_criterion": first["criterion"] if first else "",
                "first_failing_t_s": first["t_s"] if first else None,
                "first_failing_step": first["step"] if first else None,
                "failures": "; ".join(failures),
                "body_motion_recorded": (condition_dir / "trajectory.csv").is_file(),
                "flagella_motion_recorded": (
                    condition_dir / "state_archive.npz"
                ).is_file(),
            }
        )
    with (output_dir / "issue61_summary.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    all_pass = all(bool(record["strict_pass"]) for record in records)
    decision = {
        "kind": "issue61_2015_1tau_torque_stability",
        "status": "pass" if all_pass else "fail",
        "run_root": str(run_root),
        "conditions": len(records),
        "strict_pass_count": sum(bool(record["strict_pass"]) for record in records),
        "reference_evidence": manifest["reference_evidence"],
        "reference_evidence_status": (
            "verified"
            if manifest["reference_evidence"]
            else "none_available_at_run_start"
        ),
        "handoff": (
            "eligible_for_followup_evaluation_only; no supported-profile promotion or 10tau-stability claim"
            if all_pass
            else "blocked; do not promote the 2015 profile or hand off to Issue #184"
        ),
        "non_goals": [
            "10tau stability",
            "dt_star adoption",
            "dataset adoption",
            "canonical model freeze",
        ],
    }
    (output_dir / "issue61_decision.json").write_text(
        json.dumps(decision, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return output_dir


def _read_json_or_yaml(path: Path) -> dict[str, Any]:
    import yaml

    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise ValueError(f"mapping required: {path}")
    return data


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument(
        "--threshold-contract",
        type=Path,
        default=Path("conf/phase2_validation/2015_stage_a_thresholds.yaml"),
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    print(
        analyze(
            run_root=args.run_root,
            threshold_contract=args.threshold_contract,
            output_dir=args.output_dir,
        )
    )
