from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import pytest

from sim_swim.analysis.model_development_evaluation import (
    build_evaluation,
    collect_rows,
)
from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    load_yaml,
    normalize_campaign_config,
)


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "conf/phase2_multi_run/2010_hex_project_long_duration_2s_issue245.yaml"


def _summary() -> dict:
    metrics = {
        "hook_angle_err_max_deg": {"max": 61.0},
        "local_attach_first_rel_err": {"max": 0.01},
        "hook_len_rel_err_max": {"max": 0.01},
        "flag_bond_rel_err_max": {"max": 0.03},
        "flag_bend_err_max_deg": {"max": 1.0},
        "flag_torsion_err_max_deg": {"max": 1.0},
        "body_spring_max_stretch_ratio": {"max": 0.04},
        "motor_force_balance_residual_ratio": {"max": 1.0e-12},
        "motor_torque_balance_residual_ratio": {"max": 0.01},
    }
    return {
        "execution": {"status": "completed"},
        "gates": {
            "finite": {"any_fail": False},
            "shape_body": {"any_fail": False},
            "shape_nonbody": {
                "any_fail": True,
                "first_failure_category": "hook",
                "first_failure_t_s": 0.00004,
            },
        },
        "all_step_metrics": metrics,
    }


def _write_long_run(root: Path) -> Path:
    config = normalize_campaign_config(load_yaml(CONFIG))
    run_dir = root / "synchronized"
    run_dir.mkdir(parents=True)
    records, rows = [], []
    for condition in build_campaign_conditions(config):
        directory = run_dir / "conditions" / condition["condition_id"]
        directory.mkdir(parents=True)
        (directory / "run_summary.json").write_text(json.dumps(_summary()))
        (directory / "performance.json").write_text(
            json.dumps({"wall_time_s": 10.0, "steps_per_s": 100.0})
        )
        (directory / "state_archive.npz").write_bytes(b"portable archive")
        record = {
            **condition,
            "output_dir": f"/cs10/original/{condition['condition_id']}",
            "time": {"dt_star": 1.0e-3, "dt_internal_s": 4.0e-5, "duration_s": 2.0},
        }
        records.append(record)
        rows.append({"condition_id": condition["condition_id"], "completed": "True"})
    with (run_dir / "summary.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["condition_id", "completed"])
        writer.writeheader()
        writer.writerows(rows)
    (run_dir / "run_manifest.json").write_text(
        json.dumps(
            {
                "base_config": "conf/sim_swim_2010_hex.yaml",
                "campaign_config": str(CONFIG.relative_to(ROOT)),
                "model_profile": load_yaml(ROOT / "conf/sim_swim_2010_hex.yaml")[
                    "model_profile"
                ],
                "git": {"commit": "test", "is_clean": True},
                "conditions": records,
            }
        )
    )
    return run_dir


def test_issue245_static_condition_contract() -> None:
    config = normalize_campaign_config(load_yaml(CONFIG))
    conditions = build_campaign_conditions(config)
    assert len(conditions) == 18
    assert all(
        c["axis_values"]["phase_seed"] == 0
        for c in conditions
        if c["axis_values"]["n_flagella"] < 6
    )
    assert {
        (c["axis_values"]["attach_seed"], c["axis_values"]["phase_seed"])
        for c in conditions
        if c["axis_values"]["n_flagella"] == 6
    } == {(0, 0), (0, 1), (0, 2)}
    overrides = config["base_overrides"]
    assert overrides["time"]["duration"] == {"value": 2.0, "unit": "s"}
    assert overrides["time"]["integration"]["dt_star"] == 1.0e-3
    assert overrides["motor"]["torque_Nm"] == 2.5e-20
    assert overrides["output"]["archive_interval_s"] == 0.001
    assert config["output"]["save_state_archive"] is True


def test_long_duration_collects_portable_archives_and_window_qc(tmp_path: Path) -> None:
    run_dir = _write_long_run(tmp_path)
    rows, _ = collect_rows(config=load_yaml(CONFIG), run_dirs=[run_dir])
    assert len(rows) == 18
    assert {row["screen_status"] for row in rows} == {"pass"}
    assert sum(bool(row["full_ring_rotation_equivalent"]) for row in rows) == 3
    assert all(len(row["state_archive_sha256"]) == 64 for row in rows)
    outputs = build_evaluation(
        config_path=CONFIG, run_dirs=[run_dir], output_dir=tmp_path / "evaluation"
    )
    assert outputs["window_qc_csv"].is_file()
    manifest = json.loads(outputs["manifest"].read_text())
    assert manifest["full_ring_rotation_equivalent_condition_ids"] == [
        "nf06__as000__ps000",
        "nf06__as000__ps001",
        "nf06__as000__ps002",
    ]


def test_long_duration_rejects_missing_or_mismatched_archive(tmp_path: Path) -> None:
    run_dir = _write_long_run(tmp_path)
    archive = run_dir / "conditions/nf01__as000__ps000/state_archive.npz"
    archive.unlink()
    with pytest.raises(FileNotFoundError, match="state_archive.npz"):
        collect_rows(config=load_yaml(CONFIG), run_dirs=[run_dir])

    run_dir = _write_long_run(tmp_path / "again")
    manifest_path = run_dir / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    record = manifest["conditions"][0]
    archive = run_dir / "conditions" / record["condition_id"] / "state_archive.npz"
    record["artifact_sha256"] = {
        "state_archive.npz": hashlib.sha256(b"wrong").hexdigest()
    }
    manifest_path.write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        collect_rows(config=load_yaml(CONFIG), run_dirs=[run_dir])
