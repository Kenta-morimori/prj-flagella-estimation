from __future__ import annotations

import csv
import json
from pathlib import Path

from sim_swim.analysis.qualification import (
    compare_campaigns,
    compare_parallel_job,
    write_report,
)


def _campaign(
    root: Path, *, hook_error: float = 0.2, completed_steps: int = 101
) -> Path:
    root.mkdir(parents=True)
    manifest = {
        "kind": "stage_a_2015",
        "stage": "motor_off",
        "duration_tau": 0.1,
        "dt_star": 1.0e-5,
        "motor_enabled": False,
        "base_config": "conf/sim_swim_2015.yaml",
        "condition_order": ["project", "paper"],
        "git": {"commit": "abc123", "is_clean": True},
        "summary_csv": str(root / "summary.csv"),
    }
    (root / "run_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    fields = [
        "profile",
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
        "max_hook_angle_err_deg",
        "max_net_force_residual_ratio",
    ]
    rows = []
    for profile in ("project", "paper"):
        rows.append(
            {
                "profile": profile,
                "status": "completed",
                "expected_steps": "101",
                "completed_steps": str(completed_steps),
                "completion_pass": str(completed_steps == 101),
                "finite_pass_all": "True",
                "final_shape_pass_nonbody": "True",
                "body_shape_pass": "True",
                "shape_pass": "True",
                "first_fail_category_nonbody": "",
                "body_fail_category": "none",
                "max_hook_angle_err_deg": str(hook_error),
                "max_net_force_residual_ratio": "1e-10",
            }
        )
    with (root / "summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return root


def _generic_campaign(root: Path, *, completed_steps: int = 251) -> Path:
    root.mkdir(parents=True)
    condition_dir = root / "baseline"
    condition_dir.mkdir()
    run_summary = condition_dir / "run_summary.json"
    run_summary.write_text(
        json.dumps(
            {
                "execution": {
                    "status": "completed",
                    "row_count": completed_steps,
                    "expected_total_steps": 251,
                    "step_indices_contiguous_from_zero": True,
                }
            }
        ),
        encoding="utf-8",
    )
    manifest = {
        "kind": "shape_stability_grid",
        "base_config": "conf/sim_swim_2010.yaml",
        "condition_order": ["baseline"],
        "git": {"commit": "abc123", "is_clean": True},
        "args": {
            "mode": "preset",
            "duration_s": 0.001,
            "dt_star": 1.0e-4,
            "torque_nm": 2.5e-20,
            "n_flagella": 3,
            "attach_seed": 0,
            "phase_seed": 0,
        },
        "summary_csv": str(root / "summary.csv"),
        "conditions": [
            {"condition_id": "baseline", "run_summary_json": str(run_summary)}
        ],
    }
    (root / "run_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    with (root / "summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "condition_id",
                "final_shape_pass_nonbody",
                "body_shape_pass",
                "first_fail_category_nonbody",
                "max_hook_len_rel_err",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "condition_id": "baseline",
                "final_shape_pass_nonbody": "True",
                "body_shape_pass": "True",
                "first_fail_category_nonbody": "none",
                "max_hook_len_rel_err": "0.01",
            }
        )
    return root


def test_campaign_comparison_accepts_cross_platform_tolerance(tmp_path: Path) -> None:
    report = compare_campaigns(
        _campaign(tmp_path / "mac"), _campaign(tmp_path / "cs10", hook_error=0.2000001)
    )

    assert report["status"] == "PASS"
    assert all(check["status"] == "pass" for check in report["checks"])


def test_campaign_comparison_rejects_incomplete_or_provenance_mismatch(
    tmp_path: Path,
) -> None:
    left = _campaign(tmp_path / "mac")
    right = _campaign(tmp_path / "cs10", completed_steps=100)
    manifest = json.loads((right / "run_manifest.json").read_text(encoding="utf-8"))
    manifest["git"]["commit"] = "other"
    (right / "run_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")

    report = compare_campaigns(left, right)

    assert report["status"] == "FAIL"
    failed = {check["name"] for check in report["checks"] if check["status"] == "fail"}
    assert {"manifest.git.commit", "project.completed_steps"} <= failed


def test_generic_campaign_comparison_requires_completed_run_summary(
    tmp_path: Path,
) -> None:
    report = compare_campaigns(
        _generic_campaign(tmp_path / "mac"),
        _generic_campaign(tmp_path / "cs10", completed_steps=250),
    )

    assert report["status"] == "FAIL"
    failed = {check["name"] for check in report["checks"] if check["status"] == "fail"}
    assert "right.baseline.run_summary.expected_total_steps" in failed


def test_parallel_comparison_checks_launcher_and_each_child(tmp_path: Path) -> None:
    serial_a = _campaign(tmp_path / "serial_a")
    serial_b = _campaign(tmp_path / "serial_b")
    parallel_a = _campaign(tmp_path / "parallel_a")
    parallel_b = _campaign(tmp_path / "parallel_b")
    job_manifest = tmp_path / "job_manifest.json"
    job_manifest.write_text(
        json.dumps(
            {
                "status": "succeeded",
                "failed_configs": [],
                "configs": [
                    {
                        "config": "conf/phase2_sweeps/a.yaml",
                        "status": "succeeded",
                        "output_dir": str(parallel_a),
                    },
                    {
                        "config": "conf/phase2_sweeps/b.yaml",
                        "status": "succeeded",
                        "output_dir": str(parallel_b),
                    },
                ],
            }
        ),
        encoding="utf-8",
    )

    report = compare_parallel_job(
        job_manifest,
        {"conf/phase2_sweeps/a.yaml": serial_a, "conf/phase2_sweeps/b.yaml": serial_b},
    )
    json_path, csv_path = write_report(report, tmp_path / "report")

    assert report["status"] == "PASS"
    assert json_path.is_file()
    assert csv_path.is_file()
    assert (tmp_path / "report" / "manifest.json").is_file()
    assert (tmp_path / "report" / "run.log").is_file()


def test_parallel_comparison_rejects_partial_failure(tmp_path: Path) -> None:
    serial = _campaign(tmp_path / "serial")
    parallel = _campaign(tmp_path / "parallel")
    job_manifest = tmp_path / "job_manifest.json"
    job_manifest.write_text(
        json.dumps(
            {
                "status": "failed",
                "failed_configs": [1],
                "configs": [
                    {
                        "config": "conf/phase2_sweeps/a.yaml",
                        "status": "failed",
                        "output_dir": str(parallel),
                    }
                ],
            }
        ),
        encoding="utf-8",
    )

    report = compare_parallel_job(job_manifest, {"conf/phase2_sweeps/a.yaml": serial})

    assert report["status"] == "FAIL"
    assert {
        check["name"] for check in report["checks"] if check["status"] == "fail"
    } >= {
        "job.status",
        "job.failed_configs",
        "conf/phase2_sweeps/a.yaml.launcher_status",
    }
