"""Assess the clean Issue #184 2015 project nf1--6 campaign."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

from sim_swim.analysis.issue61_2015_1tau import (
    _body_drift_metrics,
    _canonical_threshold_row,
    _condition_dir,
    _first_gate_failure,
    _first_threshold_crossing,
    _float,
    _read_json,
    _read_json_or_yaml,
    _rows,
)


EXPECTED_IDS = tuple(f"nf{index:02d}" for index in range(1, 7))
EXPECTED_TORQUE_NM = 2.5e-20


def _close(actual: Any, expected: float) -> bool:
    return math.isclose(_float(actual), expected, rel_tol=1e-12, abs_tol=1e-30)


def _validate_manifest(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    if manifest.get("kind") != "generic_multi_run" or not manifest.get(
        "parallel_aggregate"
    ):
        raise ValueError("Issue #184 requires an aggregated generic campaign")
    if tuple(manifest.get("condition_order", ())) != EXPECTED_IDS:
        raise ValueError("Issue #184 requires exactly nf01 through nf06 in order")
    conditions = manifest.get("conditions")
    if not isinstance(conditions, list) or len(conditions) != len(EXPECTED_IDS):
        raise ValueError("Issue #184 requires six condition manifests")
    by_id = {str(item.get("condition_id")): item for item in conditions}
    if tuple(by_id) != EXPECTED_IDS:
        raise ValueError("Issue #184 condition manifests are incomplete or unordered")
    for index, condition_id in enumerate(EXPECTED_IDS, start=1):
        condition = by_id[condition_id]
        overrides = dict(condition.get("config_overrides", {}) or {})
        if overrides.get("flagella", {}).get("placement_mode") != "seeded_surface":
            raise ValueError(f"{condition_id} is not seeded_surface")
        if (
            overrides.get("seed", {}).get("attach_seed") != 0
            or overrides.get("seed", {}).get("phase_seed") != 0
        ):
            raise ValueError(f"{condition_id} does not use attach/phase seed 0")
        if overrides.get("time", {}).get("scale_policy") != "reference_torque":
            raise ValueError(f"{condition_id} lacks reference_torque policy")
        motor = dict(overrides.get("motor", {}) or {})
        if not all(
            _close(motor.get(name), EXPECTED_TORQUE_NM)
            for name in ("torque_Nm", "reference_torque_Nm")
        ):
            raise ValueError(f"{condition_id} torque contract mismatch")
        if condition.get("axis_values", {}).get("n_flagella") != index:
            raise ValueError(f"{condition_id} n_flagella contract mismatch")
        geometry = dict(condition.get("geometry_preflight", {}) or {})
        attachments = geometry.get("attachments")
        if (
            geometry.get("placement_mode") != "seeded_surface"
            or geometry.get("attach_seed") != 0
            or geometry.get("phase_seed") != 0
            or not isinstance(attachments, list)
            or len(attachments) != index
            or len({item.get("body_bead_index") for item in attachments}) != index
        ):
            raise ValueError(f"{condition_id} geometry preflight mismatch")
    return [by_id[condition_id] for condition_id in EXPECTED_IDS]


def analyze(*, run_root: Path, threshold_contract: Path, output_dir: Path) -> Path:
    """Write the Issue #184 diagnostic-only strict-QC decision."""
    manifest = _read_json(run_root / "run_manifest.json")
    conditions = _validate_manifest(manifest)
    completion = _read_json(run_root / "campaign_completion.json")
    if completion.get("status") != "completed" or completion.get("exit_code") != 0:
        raise ValueError("Issue #184 aggregate campaign is incomplete")
    threshold_data = _read_json_or_yaml(threshold_contract)
    thresholds = threshold_data.get("thresholds")
    if threshold_data.get("status") != "locked" or not isinstance(thresholds, dict):
        raise ValueError("threshold contract must be locked and contain thresholds")
    rows_by_id = {
        row.get("condition_id"): row for row in _rows(run_root / "summary.csv")
    }
    if set(rows_by_id) != set(EXPECTED_IDS):
        raise ValueError("Issue #184 summary.csv must contain exactly six conditions")
    output_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, Any]] = []
    for condition in conditions:
        condition_id = str(condition["condition_id"])
        condition_dir = _condition_dir(run_root, condition)
        run_summary = _read_json(condition_dir / "run_summary.json")
        gate_failure = _first_gate_failure(run_summary)
        observed = _canonical_threshold_row(rows_by_id[condition_id])
        for metric, value in _body_drift_metrics(condition_dir).items():
            if not math.isfinite(_float(observed.get(metric))):
                observed[metric] = str(value)
        failures = [
            metric
            for metric, limit in thresholds.items()
            if not math.isfinite(_float(observed.get(metric)))
            or _float(observed.get(metric)) > _float(limit)
        ]
        crossings = [
            _first_threshold_crossing(
                condition_dir, criterion=metric, limit=_float(thresholds[metric])
            )
            for metric in failures
        ]
        observed_crossings = [item for item in crossings if item is not None]
        first = gate_failure or (
            min(observed_crossings, key=lambda item: int(item["step"]))
            if observed_crossings
            else (
                {"criterion": failures[0], "t_s": None, "step": None}
                if failures
                else None
            )
        )
        row = rows_by_id[condition_id]
        records.append(
            {
                "condition_id": condition_id,
                "n_flagella": condition["axis_values"]["n_flagella"],
                "wall_time_s": _float(row.get("wall_time_s")),
                "steps_per_s": _float(row.get("steps_per_s")),
                "strict_pass": not (gate_failure or failures),
                "first_failing_criterion": first["criterion"] if first else "",
                "first_failing_t_s": first["t_s"] if first else None,
                "first_failing_step": first["step"] if first else None,
                "failures": "; ".join(
                    ([gate_failure["criterion"]] if gate_failure else []) + failures
                ),
                "body_motion_recorded": (condition_dir / "trajectory.csv").is_file(),
                "flagella_motion_recorded": (
                    condition_dir / "state_archive.npz"
                ).is_file(),
            }
        )
    with (output_dir / "issue184_summary.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    all_pass = all(bool(record["strict_pass"]) for record in records)
    (output_dir / "issue184_decision.json").write_text(
        json.dumps(
            {
                "kind": "issue184_2015_project_nf1_6_10tau",
                "status": "pass" if all_pass else "fail",
                "run_root": str(run_root),
                "conditions": len(records),
                "strict_pass_count": sum(
                    bool(record["strict_pass"]) for record in records
                ),
                "scope": "diagnostic_only",
                "handoff": "not eligible for dataset adoption, profile promotion, canonical selection, or Phase 3 handoff",
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return output_dir


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--threshold-contract",
        type=Path,
        default=Path("conf/phase2_validation/2015_stage_a_thresholds.yaml"),
    )
    args = parser.parse_args(argv)
    print(
        analyze(
            run_root=args.run_root,
            threshold_contract=args.threshold_contract,
            output_dir=args.output_dir,
        )
    )
