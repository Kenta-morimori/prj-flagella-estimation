"""Evidence-backed runtime projection for completed parallel generic campaigns."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _positive(value: Any, name: str) -> float:
    result = float(value)
    if not math.isfinite(result) or result <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return result


def estimate_parallel_runtime(
    job_root: Path,
    *,
    target_duration_s: float,
    expected_conditions: tuple[str, ...],
    historical_performance: dict[str, Path] | None = None,
) -> dict[str, Any]:
    """Project target cost from complete shards; never infer physical acceptance."""
    root = job_root.resolve()
    target = _positive(target_duration_s, "target_duration_s")
    job_path = root / "job_manifest.json"
    job = _json(job_path)
    if job.get("status") != "succeeded" or job.get("failed_configs"):
        raise ValueError("parallel job is not fully succeeded")
    aggregation = job.get("aggregation", {})
    if aggregation.get("status") != "completed":
        raise ValueError("parallel aggregate is not completed")
    campaign = root / "campaign"
    completion = _json(campaign / "campaign_completion.json")
    if completion.get("status") != "completed" or completion.get("exit_code") != 0:
        raise ValueError("campaign completion is not successful")
    manifest_path = campaign / "run_manifest.json"
    manifest = _json(manifest_path)
    records = job.get("configs")
    condition_records = manifest.get("conditions")
    if not isinstance(records, list) or not isinstance(condition_records, list):
        raise ValueError("condition records are missing")
    if [item.get("condition_id") for item in records] != list(expected_conditions):
        raise ValueError("job condition order/count mismatch")
    if [item.get("condition_id") for item in condition_records] != list(
        expected_conditions
    ):
        raise ValueError("campaign condition order/count mismatch")
    if any(item.get("status") != "succeeded" for item in records):
        raise ValueError("one or more child shards did not succeed")
    workers = int(job.get("execution", {}).get("max_workers", 0))
    if workers < 1:
        raise ValueError("invalid worker count")
    rows: list[dict[str, Any]] = []
    for index, (job_record, condition) in enumerate(
        zip(records, condition_records, strict=True), start=1
    ):
        condition_id = str(condition["condition_id"])
        expected_suffix = Path("children") / f"{index:03d}_{condition_id}" / "run"
        if not str(job_record.get("output_dir", "")).endswith(str(expected_suffix)):
            raise ValueError(f"child output path mismatch: {condition_id}")
        child_root = root / expected_suffix
        child_manifest = _json(child_root / "run_manifest.json")
        child_conditions = child_manifest.get("conditions")
        if not isinstance(child_conditions, list) or len(child_conditions) != 1:
            raise ValueError(f"child manifest mismatch: {condition_id}")
        if child_conditions[0].get("condition_id") != condition_id:
            raise ValueError(f"child condition mismatch: {condition_id}")
        time = condition.get("time")
        child_time = child_conditions[0].get("time")
        if not isinstance(time, dict) or time != child_time:
            raise ValueError(f"child time/torque contract mismatch: {condition_id}")
        duration = _positive(time.get("duration_s"), "duration_s")
        if duration >= target:
            raise ValueError("probe must be shorter than target")
        total_steps = int(time.get("total_steps", 0))
        if total_steps <= 0 or time.get("time_scale_policy") != "reference_torque":
            raise ValueError(f"invalid step/time policy: {condition_id}")
        if not math.isclose(
            float(time.get("motor_torque_Nm", 0)),
            float(time.get("reference_torque_Nm", 0)),
            rel_tol=1e-12,
        ):
            raise ValueError(f"tracking torque mismatch: {condition_id}")
        condition_dir = campaign / "conditions" / condition_id
        if (
            not condition_dir.is_symlink()
            or condition_dir.resolve() != (child_root / condition_id).resolve()
        ):
            raise ValueError(f"condition symlink mismatch: {condition_id}")
        performance_path = condition_dir / "performance.json"
        perf = _json(performance_path)
        summary_path = condition_dir / "run_summary.json"
        summary = _json(summary_path)
        completed = int(perf.get("completed_steps", -1))
        if completed != total_steps or int(perf.get("total_steps", -1)) != total_steps:
            raise ValueError(f"performance step mismatch: {condition_id}")
        if summary.get("execution", {}).get("status") != "completed":
            raise ValueError(f"condition summary is incomplete: {condition_id}")
        wall = _positive(perf.get("wall_time_s"), "wall_time_s")
        factor = target / duration
        historical = (historical_performance or {}).get(condition_id)
        historical_target = None
        historical_path = None
        historical_sha = None
        if historical is not None:
            old = _json(historical)
            old_steps = int(old.get("completed_steps", -1))
            if old_steps <= 0 or old_steps != int(old.get("total_steps", -2)):
                raise ValueError(f"historical run is incomplete: {condition_id}")
            historical_target = _positive(old.get("wall_time_s"), "historical wall") * (
                target / (duration * old_steps / total_steps)
            )
            historical_path = str(historical)
            historical_sha = _hash(historical)
        rows.append(
            {
                "condition_id": condition_id,
                "probe_duration_s": duration,
                "probe_steps": total_steps,
                "probe_wall_time_s": wall,
                "probe_steps_per_s": completed / wall,
                "target_duration_s": target,
                "target_steps": round(total_steps * factor),
                "extrapolation_factor": factor,
                "projected_wall_time_s": wall * factor,
                "historical_projected_wall_time_s": historical_target,
                "projection_to_historical_ratio": (wall * factor / historical_target)
                if historical_target
                else None,
                "strict_qc_status": "not_evaluated; runtime completion is not physical PASS",
                "online_shape_nonbody_any_fail": summary.get("gates", {})
                .get("shape_nonbody", {})
                .get("any_fail"),
                "online_shape_body_any_fail": summary.get("gates", {})
                .get("shape_body", {})
                .get("any_fail"),
                "source_git_commit": child_manifest.get("git", {}).get("commit"),
                "source_child_output_path": job_record["output_dir"],
                "performance_path": str(performance_path),
                "performance_sha256": _hash(performance_path),
                "run_summary_sha256": _hash(summary_path),
                "historical_performance_path": historical_path,
                "historical_performance_sha256": historical_sha,
            }
        )
    # Independent workers take the next condition in config order as soon as free.
    availability = [0.0] * workers
    for row in rows:
        slot = min(range(workers), key=lambda index: availability[index])
        availability[slot] += row["projected_wall_time_s"]
    return {
        "kind": "parallel_runtime_projection",
        "status": "completed",
        "decision_scope": "runtime_evidence_only",
        "job_id": job.get("job_id"),
        "job_manifest_path": str(job_path),
        "job_manifest_sha256": _hash(job_path),
        "campaign_manifest_path": str(manifest_path),
        "campaign_manifest_sha256": _hash(manifest_path),
        "target_duration_s": target,
        "worker_count": workers,
        "projected_makespan_s": max(availability),
        "projected_worker_time_s": sum(availability),
        "schedule_assumption": "fixed config order; each worker runs one condition at a time; linear step scaling; no contention correction",
        "uncertainty": "short-run startup, I/O, contention and nonlinear runtime may make the extrapolation inaccurate; historical ratios are calibration only",
        "conditions": rows,
    }


def write_runtime_projection(result: dict[str, Any], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "runtime_projection.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    rows = result["conditions"]
    with (output_dir / "runtime_projection.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
