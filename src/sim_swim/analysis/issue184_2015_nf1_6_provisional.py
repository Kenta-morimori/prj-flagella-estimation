"""Summarize explicitly-provenanced mixed-topology Issue #184 evidence."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

from sim_swim.analysis.issue184_2015_nf1_6 import EXPECTED_IDS, EXPECTED_TORQUE_NM
from sim_swim.analysis.issue61_2015_1tau import (
    _body_drift_metrics,
    _canonical_threshold_row,
    _first_gate_failure,
    _first_threshold_crossing,
    _float,
    _read_json,
    _read_json_or_yaml,
    _rows,
)


EXPECTED_PLACEMENT = {
    **{f"nf{index:02d}": "seeded_center_layer" for index in (1, 2, 3, 6)},
    "nf04": "seeded_surface",
    "nf05": "seeded_surface",
}


def _close(actual: Any, expected: float) -> bool:
    return math.isclose(_float(actual), expected, rel_tol=1e-12, abs_tol=1e-30)


def _parse_sources(values: list[str]) -> dict[str, Path]:
    sources: dict[str, Path] = {}
    for value in values:
        condition_id, separator, raw_path = value.partition("=")
        if not separator or condition_id not in EXPECTED_IDS or not raw_path:
            raise ValueError("--source must be CONDITION_ID=RUN_ROOT")
        if condition_id in sources:
            raise ValueError(f"duplicate source for {condition_id}")
        sources[condition_id] = Path(raw_path)
    if tuple(sources) != EXPECTED_IDS:
        raise ValueError(
            "provisional analysis requires exactly nf01 through nf06 in order"
        )
    return sources


def _condition_from_source(
    root: Path, condition_id: str
) -> tuple[dict[str, Any], dict[str, Any], Path]:
    manifest = _read_json(root / "run_manifest.json")
    if manifest.get("kind") != "generic_multi_run":
        raise ValueError(f"{condition_id} source is not a generic multi-run artifact")
    conditions = manifest.get("conditions")
    if not isinstance(conditions, list) or len(conditions) != 1:
        raise ValueError(f"{condition_id} source must contain exactly one condition")
    condition = conditions[0]
    if condition.get("condition_id") != condition_id:
        raise ValueError(f"source condition mismatch: expected {condition_id}")
    condition_dir = Path(str(condition.get("output_dir", "")))
    if not condition_dir.is_dir():
        raise ValueError(f"{condition_id} source condition output is missing")
    summary = _read_json(condition_dir / "run_summary.json")
    if summary.get("execution", {}).get("status") != "completed":
        raise ValueError(f"{condition_id} source is incomplete")
    return manifest, condition, condition_dir


def _validate_source(
    *,
    root: Path,
    condition_id: str,
    manifest: dict[str, Any],
    condition: dict[str, Any],
) -> dict[str, Any]:
    index = int(condition_id[-2:])
    overrides = dict(condition.get("config_overrides", {}) or {})
    placement = overrides.get("flagella", {}).get("placement_mode")
    if placement != EXPECTED_PLACEMENT[condition_id]:
        raise ValueError(f"{condition_id} placement mode mismatch: {placement}")
    seed = dict(overrides.get("seed", {}) or {})
    if seed.get("attach_seed") != 0 or seed.get("phase_seed") != 0:
        raise ValueError(f"{condition_id} does not use attach/phase seed 0")
    time = dict(overrides.get("time", {}) or {})
    if time.get("scale_policy") != "reference_torque":
        raise ValueError(f"{condition_id} lacks reference_torque policy")
    if not _close(time.get("duration", {}).get("value"), 10.0):
        raise ValueError(f"{condition_id} duration contract mismatch")
    if not _close(time.get("integration", {}).get("dt_star"), 1.0e-5):
        raise ValueError(f"{condition_id} dt_star contract mismatch")
    motor = dict(overrides.get("motor", {}) or {})
    if not all(
        _close(motor.get(name), EXPECTED_TORQUE_NM)
        for name in ("torque_Nm", "reference_torque_Nm")
    ):
        raise ValueError(f"{condition_id} torque contract mismatch")
    if condition.get("axis_values", {}).get("n_flagella") != index:
        raise ValueError(f"{condition_id} n_flagella contract mismatch")
    topology = (
        dict(condition.get("geometry", {}) or {})
        .get("actual", {})
        .get("attachment_topology")
    )
    if (
        not isinstance(topology, list)
        or len(topology) != index
        or len({item.get("body_bead_index") for item in topology}) != index
    ):
        raise ValueError(f"{condition_id} attachment topology mismatch")
    return {
        "condition_id": condition_id,
        "source_root": str(root),
        "git_commit": manifest.get("git", {}).get("commit"),
        "placement_mode": placement,
        "attach_seed": seed["attach_seed"],
        "phase_seed": seed["phase_seed"],
        "attachment_topology": topology,
    }


def analyze(
    *, sources: dict[str, Path], threshold_contract: Path, output_dir: Path
) -> Path:
    """Write a diagnostic-only report for six explicitly selected child runs."""
    if tuple(sources) != EXPECTED_IDS:
        raise ValueError(
            "provisional analysis requires exactly nf01 through nf06 in order"
        )
    threshold_data = _read_json_or_yaml(threshold_contract)
    thresholds = threshold_data.get("thresholds")
    if threshold_data.get("status") != "locked" or not isinstance(thresholds, dict):
        raise ValueError("threshold contract must be locked and contain thresholds")
    output_dir.mkdir(parents=True, exist_ok=True)
    records: list[dict[str, Any]] = []
    provenance: list[dict[str, Any]] = []
    for condition_id in EXPECTED_IDS:
        root = sources[condition_id]
        manifest, condition, condition_dir = _condition_from_source(root, condition_id)
        source = _validate_source(
            root=root, condition_id=condition_id, manifest=manifest, condition=condition
        )
        provenance.append(source)
        row_by_id = {
            row.get("condition_id"): row for row in _rows(root / "summary.csv")
        }
        if set(row_by_id) != {condition_id}:
            raise ValueError(f"{condition_id} source summary must contain only itself")
        observed = _canonical_threshold_row(row_by_id[condition_id])
        for metric, value in _body_drift_metrics(condition_dir).items():
            if not math.isfinite(_float(observed.get(metric))):
                observed[metric] = str(value)
        run_summary = _read_json(condition_dir / "run_summary.json")
        gate_failure = _first_gate_failure(run_summary)
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
        records.append(
            {
                "condition_id": condition_id,
                "n_flagella": int(condition_id[-2:]),
                "source_root": str(root),
                "git_commit": source["git_commit"],
                "placement_mode": source["placement_mode"],
                "attachment_topology": json.dumps(source["attachment_topology"]),
                "wall_time_s": _float(row_by_id[condition_id].get("wall_time_s")),
                "steps_per_s": _float(row_by_id[condition_id].get("steps_per_s")),
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
    with (output_dir / "issue184_provisional_summary.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    (output_dir / "provenance_manifest.json").write_text(
        json.dumps(
            {
                "kind": "issue184_2015_project_nf1_6_provisional_provenance",
                "scope": "mixed_topology_diagnostic_only",
                "conditions": provenance,
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    strict_pass_count = sum(bool(record["strict_pass"]) for record in records)
    (output_dir / "issue184_provisional_decision.json").write_text(
        json.dumps(
            {
                "kind": "issue184_2015_project_nf1_6_provisional",
                "status": "provisional",
                "strict_status": "pass"
                if strict_pass_count == len(records)
                else "fail",
                "strict_pass_count": strict_pass_count,
                "conditions": len(records),
                "scope": "mixed_topology_diagnostic_only",
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
    parser.add_argument(
        "--source", action="append", default=[], metavar="CONDITION_ID=RUN_ROOT"
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--threshold-contract",
        type=Path,
        default=Path("conf/phase2_validation/2015_stage_a_thresholds.yaml"),
    )
    args = parser.parse_args(argv)
    print(
        analyze(
            sources=_parse_sources(args.source),
            threshold_contract=args.threshold_contract,
            output_dir=args.output_dir,
        )
    )
