"""Completed-shard and scheduling checks for runtime cost evidence."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from sim_swim.analysis.runtime_projection import estimate_parallel_runtime


def _write(path: Path, value: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value), encoding="utf-8")


def _campaign(tmp_path: Path) -> Path:
    root = tmp_path / "job"
    campaign = root / "campaign"
    conditions = ["nf01", "nf02", "nf03"]
    job_records = []
    campaign_records = []
    for index, name in enumerate(conditions, start=1):
        child = root / "children" / f"{index:03d}_{name}" / "run"
        condition_dir = child / name
        condition_dir.mkdir(parents=True)
        link = campaign / "conditions" / name
        link.parent.mkdir(parents=True, exist_ok=True)
        link.symlink_to(condition_dir)
        time = {
            "duration_s": 0.01,
            "total_steps": 25000,
            "time_scale_policy": "reference_torque",
            "motor_torque_Nm": 2.5e-20,
            "reference_torque_Nm": 2.5e-20,
        }
        _write(
            child / "run_manifest.json",
            {
                "git": {"commit": "abc"},
                "conditions": [{"condition_id": name, "time": time}],
            },
        )
        _write(
            condition_dir / "performance.json",
            {
                "wall_time_s": 100.0 * index,
                "completed_steps": 25000,
                "total_steps": 25000,
            },
        )
        _write(
            condition_dir / "run_summary.json", {"execution": {"status": "completed"}}
        )
        job_records.append(
            {"condition_id": name, "output_dir": str(child), "status": "succeeded"}
        )
        campaign_records.append({"condition_id": name, "time": time})
    _write(
        root / "job_manifest.json",
        {
            "status": "succeeded",
            "failed_configs": [],
            "aggregation": {"status": "completed"},
            "execution": {"max_workers": 2},
            "configs": job_records,
        },
    )
    _write(campaign / "run_manifest.json", {"conditions": campaign_records})
    _write(
        campaign / "campaign_completion.json", {"status": "completed", "exit_code": 0}
    )
    return root


def test_projection_uses_each_completed_shard_and_worker_schedule(
    tmp_path: Path,
) -> None:
    root = _campaign(tmp_path)
    historical = tmp_path / "old-nf01-performance.json"
    _write(
        historical,
        {"wall_time_s": 200000, "completed_steps": 1000000, "total_steps": 1000000},
    )
    result = estimate_parallel_runtime(
        root,
        target_duration_s=0.5,
        expected_conditions=("nf01", "nf02", "nf03"),
        historical_performance={"nf01": historical},
    )
    assert [row["target_steps"] for row in result["conditions"]] == [1250000] * 3
    assert [row["projected_wall_time_s"] for row in result["conditions"]] == [
        5000,
        10000,
        15000,
    ]
    assert result["projected_makespan_s"] == 20000
    assert result["projected_worker_time_s"] == 30000
    assert result["decision_scope"] == "runtime_evidence_only"
    assert result["conditions"][0]["historical_projected_wall_time_s"] == 250000
    assert result["conditions"][0]["projection_to_historical_ratio"] == 0.02


@pytest.mark.parametrize(
    "mutation",
    ["job_failed", "missing_shard", "partial_steps", "time_mismatch", "bad_link"],
)
def test_projection_rejects_incomplete_or_inconsistent_evidence(
    tmp_path: Path, mutation: str
) -> None:
    root = _campaign(tmp_path)
    job_path = root / "job_manifest.json"
    job = json.loads(job_path.read_text())
    if mutation == "job_failed":
        job["status"] = "failed"
    elif mutation == "missing_shard":
        job["configs"].pop()
    elif mutation == "partial_steps":
        perf = root / "campaign" / "conditions" / "nf02" / "performance.json"
        data = json.loads(perf.read_text())
        data["completed_steps"] = 24999
        _write(perf, data)
    elif mutation == "time_mismatch":
        manifest = root / "campaign" / "run_manifest.json"
        data = json.loads(manifest.read_text())
        data["conditions"][0]["time"]["motor_torque_Nm"] = 1e-21
        _write(manifest, data)
    elif mutation == "bad_link":
        link = root / "campaign" / "conditions" / "nf03"
        link.unlink()
        link.symlink_to(root / "children" / "001_nf01" / "run" / "nf01")
    _write(job_path, job)
    with pytest.raises(ValueError):
        estimate_parallel_runtime(
            root, target_duration_s=0.5, expected_conditions=("nf01", "nf02", "nf03")
        )
