from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]


def _load_script(name: str, rel_path: str):
    path = ROOT / rel_path
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


run_sweep = _load_script(
    "run_flagella_count_behavior_sweep",
    "scripts/analysis/run_flagella_count_behavior_sweep.py",
)
build_dataset = _load_script(
    "build_flagella_count_behavior_dataset",
    "scripts/analysis/build_flagella_count_behavior_dataset.py",
)


def test_flagella_count_conditions_use_expected_sample_ids() -> None:
    config = {
        "sweep": {
            "n_flagella": [1, 2, 3, 6],
            "seeds": [0, 1, 2],
        }
    }

    conditions = run_sweep.build_conditions(config)

    assert len(conditions) == 12
    assert conditions[0] == {
        "sample_id": "nf01_seed000",
        "condition_tag": "n_flagella=1,seed=0",
        "n_flagella": 1,
        "seed": 0,
    }
    assert conditions[-1]["sample_id"] == "nf06_seed002"


def test_flagella_count_dry_run_writes_manifest_and_configs(tmp_path: Path) -> None:
    analysis_config = {
        "dataset_id": "test_dataset",
        "run_batch_id": "test_dataset",
        "base_config": str(ROOT / "conf/sim_swim.yaml"),
        "feature_schema": str(
            ROOT / "conf/analysis/flagella_count_behavior_features.yaml"
        ),
        "sweep": {"n_flagella": [1, 2], "seeds": [0, 1]},
        "base_overrides": {
            "time.duration_s": 0.5,
            "time.dt_star": 1.0e-4,
            "motor.torque_Nm": 2.5e-20,
            "motor.enable_switching": False,
            "motor.force_distribution": "material_twist_local_couple",
            "flagella.initial_helix_axis_from_rear_deg": 0,
        },
        "output": {
            "run_batch_dir": str(tmp_path / "runs/test_dataset"),
            "dataset_dir": str(tmp_path / "datasets/test_dataset"),
            "save_numeric_timeseries": True,
            "timeseries_sampling": "all_steps",
        },
    }
    config_path = tmp_path / "analysis.yaml"
    config_path.write_text(yaml.safe_dump(analysis_config), encoding="utf-8")

    manifest_path = run_sweep.run_batch(
        analysis_config_path=config_path,
        dry_run=True,
        overwrite=False,
        stop_on_shape_fail=False,
        sample_limit=None,
        progress_interval=None,
    )

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["run_batch_id"] == "test_dataset"
    assert len(manifest["samples"]) == 4
    assert {sample["status"] for sample in manifest["samples"]} == {"planned"}
    sample_config = yaml.safe_load(
        (tmp_path / "runs/test_dataset/configs/nf02_seed001.yaml").read_text(
            encoding="utf-8"
        )
    )
    assert sample_config["flagella"]["n_flagella"] == 2
    assert sample_config["seed"]["global_seed"] == 1
    assert sample_config["time"]["dt_star"] == 1.0e-4


def _write_step_summary(path: Path, rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "step",
        "t_s",
        "dt_s",
        "finite_pass",
        "shape_pass_nonbody_strict",
        "shape_pass_nonbody_hook_len_relaxed",
        "first_fail_category_nonbody",
        "body_displacement_um",
        "body_speed_um_s",
        "body_axis_step_angle_deg",
        "body_axis_cumulative_angle_deg",
        "body_axis_wobble_rms_deg",
        "body_angular_velocity_rms_rad_s",
        "flag_helix_axis_alignment_order",
        "flag_helix_axis_mean_deviation_deg_max",
        "flag_helix_axis_pair_angle_deg_mean",
        "flag_helix_axis_pair_angle_deg_max",
        "flag_helix_axis_rearward_projection_min",
        "bundle_axis_vs_body_axis_angle_deg",
        "hook_len_rel_err_max",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def test_dataset_builder_outputs_summary_qc_and_timeseries(tmp_path: Path) -> None:
    run_dir = tmp_path / "runs/test_dataset"
    dataset_dir = tmp_path / "datasets/test_dataset"
    raw_a = run_dir / "samples/nf01_seed000/raw"
    raw_b = run_dir / "samples/nf02_seed000/raw"
    _write_step_summary(
        raw_a / "step_summary.csv",
        [
            {
                "step": 0,
                "t_s": 0.0,
                "dt_s": 0.1,
                "finite_pass": True,
                "shape_pass_nonbody_strict": True,
                "shape_pass_nonbody_hook_len_relaxed": True,
                "first_fail_category_nonbody": "none",
                "body_displacement_um": 0.0,
                "body_speed_um_s": 0.0,
                "body_axis_step_angle_deg": 0.0,
                "body_axis_cumulative_angle_deg": 0.0,
                "body_axis_wobble_rms_deg": 0.0,
                "body_angular_velocity_rms_rad_s": 0.0,
                "flag_helix_axis_alignment_order": 1.0,
                "flag_helix_axis_mean_deviation_deg_max": 0.0,
                "flag_helix_axis_pair_angle_deg_mean": 0.0,
                "flag_helix_axis_pair_angle_deg_max": 0.0,
                "flag_helix_axis_rearward_projection_min": 1.0,
                "bundle_axis_vs_body_axis_angle_deg": 5.0,
                "hook_len_rel_err_max": 0.1,
            },
            {
                "step": 1,
                "t_s": 0.1,
                "dt_s": 0.1,
                "finite_pass": True,
                "shape_pass_nonbody_strict": True,
                "shape_pass_nonbody_hook_len_relaxed": True,
                "first_fail_category_nonbody": "none",
                "body_displacement_um": 1.0,
                "body_speed_um_s": 10.0,
                "body_axis_step_angle_deg": 2.0,
                "body_axis_cumulative_angle_deg": 2.0,
                "body_axis_wobble_rms_deg": 1.0,
                "body_angular_velocity_rms_rad_s": 0.3,
                "flag_helix_axis_alignment_order": 0.95,
                "flag_helix_axis_mean_deviation_deg_max": 4.0,
                "flag_helix_axis_pair_angle_deg_mean": 5.0,
                "flag_helix_axis_pair_angle_deg_max": 8.0,
                "flag_helix_axis_rearward_projection_min": 0.98,
                "bundle_axis_vs_body_axis_angle_deg": 6.0,
                "hook_len_rel_err_max": 0.2,
            },
        ],
    )
    _write_step_summary(
        raw_b / "step_summary.csv",
        [
            {
                "step": 0,
                "t_s": 0.0,
                "dt_s": 0.1,
                "finite_pass": True,
                "shape_pass_nonbody_strict": False,
                "shape_pass_nonbody_hook_len_relaxed": True,
                "first_fail_category_nonbody": "hook",
                "body_displacement_um": 0.0,
                "body_speed_um_s": "",
                "body_axis_step_angle_deg": 0.0,
                "body_axis_cumulative_angle_deg": 0.0,
                "body_axis_wobble_rms_deg": 0.0,
                "body_angular_velocity_rms_rad_s": 0.0,
                "flag_helix_axis_alignment_order": 0.8,
                "flag_helix_axis_mean_deviation_deg_max": 12.0,
                "flag_helix_axis_pair_angle_deg_mean": 14.0,
                "flag_helix_axis_pair_angle_deg_max": 20.0,
                "flag_helix_axis_rearward_projection_min": 0.9,
                "bundle_axis_vs_body_axis_angle_deg": "",
                "hook_len_rel_err_max": 1.5,
            }
        ],
    )

    feature_schema = tmp_path / "features.yaml"
    feature_schema.write_text("feature_categories: {}\n", encoding="utf-8")
    analysis_config = {
        "dataset_id": "test_dataset",
        "feature_schema": str(feature_schema),
        "output": {
            "run_batch_dir": str(run_dir),
            "dataset_dir": str(dataset_dir),
            "timeseries_sampling": "all_steps",
        },
    }
    analysis_config_path = tmp_path / "analysis.yaml"
    analysis_config_path.write_text(yaml.safe_dump(analysis_config), encoding="utf-8")
    run_manifest = {
        "run_batch_id": "test_dataset",
        "dataset_id": "test_dataset",
        "feature_schema": str(feature_schema),
        "samples": [
            {
                "sample_id": "nf01_seed000",
                "n_flagella": 1,
                "seed": 0,
                "duration_s": 0.1,
                "dt_star": 1.0e-4,
                "torque_Nm": 2.5e-20,
                "force_distribution": "material_twist_local_couple",
                "condition_tag": "n_flagella=1,seed=0",
                "status": "completed",
                "raw_dir": str(raw_a),
            },
            {
                "sample_id": "nf02_seed000",
                "n_flagella": 2,
                "seed": 0,
                "duration_s": 0.1,
                "dt_star": 1.0e-4,
                "torque_Nm": 2.5e-20,
                "force_distribution": "material_twist_local_couple",
                "condition_tag": "n_flagella=2,seed=0",
                "status": "completed",
                "raw_dir": str(raw_b),
            },
        ],
    }
    run_manifest_path = run_dir / "run_manifest.json"
    run_manifest_path.parent.mkdir(parents=True, exist_ok=True)
    run_manifest_path.write_text(
        json.dumps(run_manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    out_dir = build_dataset.build_dataset(
        analysis_config_path=analysis_config_path,
        run_manifest_path=run_manifest_path,
        overwrite=True,
    )

    assert out_dir == dataset_dir
    assert (dataset_dir / "dataset_manifest.json").is_file()
    assert (dataset_dir / "feature_schema_used.yaml").read_text(
        encoding="utf-8"
    ) == "feature_categories: {}\n"
    with (dataset_dir / "summary.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2
    assert rows[0]["quality_class"] == "strict_pass"
    assert float(rows[0]["cell_displacement"]) == 1.0
    assert float(rows[0]["cell_straightness"]) == 1.0
    assert rows[1]["quality_class"] == "relaxed_pass"
    assert rows[1]["hook_wrapped"] == "True"
    assert (dataset_dir / "qc_summary.csv").is_file()
    with (dataset_dir / "timeseries/nf01_seed000.csv").open(
        encoding="utf-8", newline=""
    ) as handle:
        ts_rows = list(csv.DictReader(handle))
    assert len(ts_rows) == 2
    assert ts_rows[0]["dataset_id"] == "test_dataset"
