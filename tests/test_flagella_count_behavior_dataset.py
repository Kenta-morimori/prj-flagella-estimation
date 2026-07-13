from __future__ import annotations

import csv
import importlib.util
import json
import math
from pathlib import Path
import sys

import numpy as np
import pytest
import yaml

from sim_swim.analysis.flagella_count_behavior import (
    load_state_archive,
    save_state_archive,
)
from sim_swim.analysis.multi_run_campaign import (
    apply_campaign_cli_overrides,
    build_campaign_conditions,
    load_yaml,
)
from sim_swim.render.video_writer import VideoRenderResult
from sim_swim.sim.core import SimulationState
from sim_swim.sim.params import SimulationConfig


ROOT = Path(__file__).resolve().parents[1]


def _load_script(name: str, rel_path: str):
    path = ROOT / rel_path
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


build_dataset = _load_script(
    "build_dataset",
    "scripts/02_phase2_analysis/build_dataset.py",
)
render_sample = _load_script(
    "render_sample",
    "scripts/02_phase2_analysis/render_sample.py",
)
plot_distributions = _load_script(
    "plot_distributions",
    "scripts/02_phase2_analysis/plot_distributions.py",
)


def test_distribution_plot_uses_integer_n_flagella_ticks() -> None:
    base = plot_distributions.pd.DataFrame(
        {
            "n_flagella": [1, 2, 2, 3, 6],
        }
    )

    ticks = plot_distributions._n_flagella_ticks(base)

    assert ticks == [1.0, 2.0, 3.0, 6.0]
    assert plot_distributions._n_flagella_tick_labels(ticks) == ["1", "2", "3", "6"]


def test_flagella_count_diagnostic_multi_run_config_generates_expected_conditions() -> (
    None
):
    config_path = ROOT / "conf/phase2_multi_run/flagella_count_behavior_diagnostic.yaml"
    config = apply_campaign_cli_overrides(load_yaml(config_path), [])

    conditions = build_campaign_conditions(config)

    assert config["dataset"]["dataset_id"] == ("fc_nf1_2_3_6_as3_ps3_torque2p0_dur1p0")
    assert config["base_overrides"]["motor"]["torque_Nm"] == pytest.approx(2.0e-20)
    assert config["base_overrides"]["motor"]["force_distribution"] == (
        "root_torque_segment_couples"
    )
    assert len(conditions) == 36
    assert conditions[0]["condition_id"] == "as000__ps000__nf01"
    assert conditions[0]["axis_values"] == {
        "attach_seed": 0,
        "phase_seed": 0,
        "n_flagella": 1,
    }
    assert conditions[-1]["condition_id"] == "as002__ps002__nf06"


def test_flagella_count_dataset_overrides_live_under_campaign_dataset() -> None:
    config_path = ROOT / "conf/phase2_multi_run/flagella_count_behavior_diagnostic.yaml"
    config = apply_campaign_cli_overrides(
        load_yaml(config_path),
        [
            "dataset.dataset_id=cli_dataset",
            "dataset.output_dir=/tmp/cli_dataset",
            "time.duration_s=0.25",
        ],
    )
    conditions = build_campaign_conditions(config)

    assert config["dataset"]["dataset_id"] == "cli_dataset"
    assert config["dataset"]["output_dir"] == "/tmp/cli_dataset"
    assert config["base_overrides"]["time"]["duration_s"] == pytest.approx(0.25)
    assert conditions[0]["config_overrides"]["time"]["duration_s"] == pytest.approx(
        0.25
    )


def test_state_archive_round_trip(tmp_path: Path) -> None:
    states = [
        SimulationState(
            t=0.0,
            position_um=(0.0, 1.0, 2.0),
            quaternion=(0.0, 0.0, 0.0, 1.0),
            velocity_um_s=(1.0, 2.0, 3.0),
            omega_rad_s=(0.1, 0.2, 0.3),
            bead_positions_um=np.array([[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]]),
            flag_states=(0, 1),
            reverse_flagella=(1,),
        ),
        SimulationState(
            t=0.1,
            position_um=(3.0, 4.0, 5.0),
            quaternion=(0.1, 0.2, 0.3, 0.4),
            velocity_um_s=(4.0, 5.0, 6.0),
            omega_rad_s=(0.4, 0.5, 0.6),
            bead_positions_um=np.array([[2.0, 2.0, 2.0], [3.0, 3.0, 3.0]]),
            flag_states=(2, 3),
            reverse_flagella=(1,),
        ),
    ]
    archive_path = tmp_path / "state_archive.npz"

    save_state_archive(archive_path, states)
    loaded = load_state_archive(archive_path)

    assert len(loaded) == 2
    assert loaded[0].position_um == (0.0, 1.0, 2.0)
    assert loaded[1].quaternion == (0.1, 0.2, 0.3, 0.4)
    assert loaded[0].flag_states == (0, 1)
    assert loaded[1].reverse_flagella == (1,)
    assert loaded[1].bead_positions_um.shape == (2, 3)


def test_render_sample_from_archive(tmp_path: Path) -> None:
    raw_cfg = yaml.safe_load((ROOT / "conf/sim_swim.yaml").read_text(encoding="utf-8"))
    raw_cfg["flagella"]["n_flagella"] = 1
    raw_cfg["output_sampling"]["out_all_steps_3d"] = False
    raw_cfg["output_sampling"]["fps_out_3d"] = 2.0
    raw_cfg["output_sampling"]["fps_out_2d"] = 2.0
    raw_cfg["render"]["save_frames_3d"] = False
    raw_cfg["render"]["save_frames_2d"] = False
    raw_cfg["render"]["render_flagella_2d"] = False
    config_path = tmp_path / "configs/nf01_seed000.yaml"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(yaml.safe_dump(raw_cfg), encoding="utf-8")
    cfg = SimulationConfig.from_dict(raw_cfg)
    sim = render_sample.Simulator(cfg)
    sample_dir = tmp_path / "runs/test_dataset/samples/nf01_seed000"
    archive_path = sample_dir / "raw/state_archive.npz"
    state = SimulationState(
        t=0.0,
        position_um=(0.0, 0.0, 0.0),
        quaternion=(0.0, 0.0, 0.0, 1.0),
        velocity_um_s=(0.0, 0.0, 0.0),
        omega_rad_s=(0.0, 0.0, 0.0),
        bead_positions_um=sim.model.positions_m.copy() * 1.0e6,
        flag_states=(),
        reverse_flagella=(),
    )
    save_state_archive(archive_path, [state])

    out_dir = render_sample.render_sample(
        sample_dir=sample_dir,
        output_dir=tmp_path / "replays/nf01_seed000",
        config_path=config_path,
        archive_path=archive_path,
    )

    assert (out_dir / "manifest.json").is_file()
    assert (out_dir / "trajectory.csv").is_file()
    assert (out_dir / "render/swim3d.mp4").is_file()
    assert (out_dir / "render2d/projection.mp4").is_file()


def test_render_sample_defaults_to_lightweight_3d_replay(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    raw_cfg = yaml.safe_load((ROOT / "conf/sim_swim.yaml").read_text(encoding="utf-8"))
    raw_cfg["flagella"]["n_flagella"] = 1
    raw_cfg["output_sampling"]["out_all_steps_3d"] = True
    raw_cfg["output_sampling"]["fps_out_3d"] = 99.0
    raw_cfg["output_sampling"]["fps_out_2d"] = 88.0
    config_path = tmp_path / "configs/nf01_seed000.yaml"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(yaml.safe_dump(raw_cfg), encoding="utf-8")
    cfg = SimulationConfig.from_dict(raw_cfg)
    sim = render_sample.Simulator(cfg)
    sample_dir = tmp_path / "runs/test_dataset/samples/nf01_seed000"
    archive_path = sample_dir / "raw/state_archive.npz"
    save_state_archive(
        archive_path,
        [
            SimulationState(
                t=0.0,
                position_um=(0.0, 0.0, 0.0),
                quaternion=(0.0, 0.0, 0.0, 1.0),
                velocity_um_s=(0.0, 0.0, 0.0),
                omega_rad_s=(0.0, 0.0, 0.0),
                bead_positions_um=sim.model.positions_m.copy() * 1.0e6,
                flag_states=(),
                reverse_flagella=(),
            )
        ],
    )
    captured: dict[str, SimulationConfig] = {}

    def fake_save_swim_movie(states, cfg, rig, out_dir) -> VideoRenderResult:
        captured["cfg_3d"] = cfg
        out_dir.mkdir(parents=True, exist_ok=True)
        return VideoRenderResult(
            path=str(out_dir / "swim3d.mp4"),
            selected_codec="avc1",
            attempted_codecs=("avc1", "H264", "mp4v"),
            fps=99.0,
            frame_size=(1000, 1000),
            frame_count=1,
        )

    def fake_project_states(states, cfg, rig, out_dir) -> VideoRenderResult:
        captured["cfg_2d"] = cfg
        out_dir.mkdir(parents=True, exist_ok=True)
        return VideoRenderResult(
            path=str(out_dir / "projection.mp4"),
            selected_codec="mp4v",
            attempted_codecs=("avc1", "H264", "mp4v"),
            fps=88.0,
            frame_size=(256, 256),
            frame_count=1,
        )

    monkeypatch.setattr(render_sample, "save_swim_movie", fake_save_swim_movie)
    monkeypatch.setattr(render_sample, "project_states", fake_project_states)

    out_dir = render_sample.render_sample(
        sample_dir=sample_dir,
        output_dir=tmp_path / "replays/nf01_seed000",
        config_path=config_path,
        archive_path=archive_path,
    )

    assert captured["cfg_3d"].output_sampling.out_all_steps_3d is False
    assert captured["cfg_3d"].output_sampling.fps_out_3d == pytest.approx(99.0)
    assert captured["cfg_2d"].output_sampling.fps_out_2d == pytest.approx(88.0)
    manifest = json.loads((out_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["render_sampling"] == {
        "out_all_steps_3d": False,
        "fps_out_3d": 99.0,
        "fps_out_2d": 88.0,
    }
    assert manifest["render_sampling_overrides"] == {"out_all_steps_3d": False}
    assert manifest["render_video"]["render3d"]["selected_codec"] == "avc1"
    assert manifest["render_video"]["render2d"]["selected_codec"] == "mp4v"


def test_render_sample_cli_passes_sampling_options(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    captured: dict[str, object] = {}
    output_dir = tmp_path / "replay"

    def fake_render_sample(**kwargs):
        captured.update(kwargs)
        return output_dir

    monkeypatch.setattr(render_sample, "render_sample", fake_render_sample)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "render_sample.py",
            "--sample-dir",
            str(tmp_path / "sample"),
            "--output-dir",
            str(output_dir),
            "--fps-out-3d",
            "12.5",
            "--fps-out-2d",
            "8",
            "--out-all-steps-3d",
        ],
    )

    render_sample.main()

    assert captured["output_dir"] == output_dir
    assert captured["out_all_steps_3d"] is True
    assert captured["fps_out_3d"] == pytest.approx(12.5)
    assert captured["fps_out_2d"] == pytest.approx(8.0)


def test_render_dataset_renders_all_samples_into_dataset_replays(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    dataset_dir = tmp_path / "dataset"
    run_dir = tmp_path / "runs/test_dataset"
    run_manifest_path = run_dir / "run_manifest.json"
    samples = [
        {
            "sample_id": "nf01_seed000",
            "sample_dir": str(run_dir / "samples/nf01_seed000"),
            "config_path": str(run_dir / "configs/nf01_seed000.yaml"),
            "outputs": {
                "state_archive_npz": str(
                    run_dir / "samples/nf01_seed000/raw/state_archive.npz"
                )
            },
        },
        {
            "sample_id": "nf02_seed000",
            "sample_dir": str(run_dir / "samples/nf02_seed000"),
            "config_path": str(run_dir / "configs/nf02_seed000.yaml"),
            "outputs": {
                "state_archive_npz": str(
                    run_dir / "samples/nf02_seed000/raw/state_archive.npz"
                )
            },
        },
    ]
    run_manifest_path.parent.mkdir(parents=True, exist_ok=True)
    run_manifest_path.write_text(json.dumps({"samples": samples}), encoding="utf-8")
    dataset_dir.mkdir(parents=True, exist_ok=True)
    (dataset_dir / "dataset_manifest.json").write_text(
        json.dumps(
            {
                "dataset_id": "test_dataset",
                "run_batch_id": "test_dataset",
                "run_manifest": str(run_manifest_path),
            }
        ),
        encoding="utf-8",
    )
    stale = dataset_dir / "replays/nf01_seed000/stale.txt"
    stale.parent.mkdir(parents=True, exist_ok=True)
    stale.write_text("old", encoding="utf-8")
    captured: list[dict[str, object]] = []

    def fake_render_sample(**kwargs):
        assert not (kwargs["output_dir"] / "stale.txt").exists()
        captured.append(kwargs)
        kwargs["output_dir"].mkdir(parents=True, exist_ok=True)
        (kwargs["output_dir"] / "manifest.json").write_text("{}", encoding="utf-8")
        return kwargs["output_dir"]

    monkeypatch.setattr(render_sample, "render_sample", fake_render_sample)

    replay_root = render_sample.render_dataset(
        dataset_dir=dataset_dir,
        output_dir=None,
        out_all_steps_3d=False,
        fps_out_3d=12.5,
        fps_out_2d=8.0,
    )

    assert replay_root == (dataset_dir / "replays").resolve()
    assert [
        item["sample_id"]
        for item in json.loads(
            (replay_root / "manifest.json").read_text(encoding="utf-8")
        )["samples"]
    ] == ["nf01_seed000", "nf02_seed000"]
    assert len(captured) == 2
    assert captured[0]["output_dir"] == replay_root / "nf01_seed000"
    assert captured[0]["fps_out_3d"] == pytest.approx(12.5)
    assert captured[0]["fps_out_2d"] == pytest.approx(8.0)


def test_render_dataset_rejects_sample_id_that_escapes_replay_root(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    dataset_dir = tmp_path / "dataset"
    run_dir = tmp_path / "runs/test_dataset"
    run_manifest_path = run_dir / "run_manifest.json"
    outside_dir = tmp_path / "outside"
    sentinel = outside_dir / "sentinel.txt"
    outside_dir.mkdir(parents=True)
    sentinel.write_text("keep", encoding="utf-8")
    sample_dir = run_dir / "samples/bad"
    run_manifest_path.parent.mkdir(parents=True, exist_ok=True)
    run_manifest_path.write_text(
        json.dumps(
            {
                "samples": [
                    {
                        "sample_id": "../outside",
                        "sample_dir": str(sample_dir),
                        "config_path": str(run_dir / "configs/bad.yaml"),
                        "outputs": {
                            "state_archive_npz": str(
                                sample_dir / "raw/state_archive.npz"
                            )
                        },
                    }
                ]
            }
        ),
        encoding="utf-8",
    )
    dataset_dir.mkdir(parents=True, exist_ok=True)
    (dataset_dir / "dataset_manifest.json").write_text(
        json.dumps({"run_manifest": str(run_manifest_path)}),
        encoding="utf-8",
    )

    def fail_render_sample(**kwargs):
        raise AssertionError("render_sample must not be called for unsafe sample_id")

    monkeypatch.setattr(render_sample, "render_sample", fail_render_sample)

    with pytest.raises(ValueError, match="plain directory name"):
        render_sample.render_dataset(
            dataset_dir=dataset_dir,
            output_dir=tmp_path / "replays",
        )

    assert sentinel.read_text(encoding="utf-8") == "keep"


def test_render_dataset_rejects_absolute_sample_id(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="plain directory name"):
        render_sample._resolve_sample_output_dir(tmp_path / "replays", "/tmp/bad")


def test_render_dataset_rejects_backslash_sample_id(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="plain directory name"):
        render_sample._resolve_sample_output_dir(tmp_path / "replays", "bad\\name")


def test_render_sample_cli_passes_dataset_dir_options(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    captured: dict[str, object] = {}
    output_dir = tmp_path / "replays"

    def fake_render_dataset(**kwargs):
        captured.update(kwargs)
        return output_dir

    monkeypatch.setattr(render_sample, "render_dataset", fake_render_dataset)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "render_sample.py",
            "--dataset-dir",
            str(tmp_path / "dataset"),
            "--output-dir",
            str(output_dir),
            "--fps-out-3d",
            "12.5",
            "--fps-out-2d",
            "8",
            "--out-all-steps-3d",
        ],
    )

    render_sample.main()

    assert captured["dataset_dir"] == tmp_path / "dataset"
    assert captured["output_dir"] == output_dir
    assert captured["out_all_steps_3d"] is True
    assert captured["fps_out_3d"] == pytest.approx(12.5)
    assert captured["fps_out_2d"] == pytest.approx(8.0)


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
        "bundle_axis_vs_rear_angle_deg",
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
    raw_a = run_dir / "as000__ps000__nf01"
    raw_b = run_dir / "as000__ps000__nf02"
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
                "bundle_axis_vs_rear_angle_deg": 6.0,
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
                "bundle_axis_vs_rear_angle_deg": 7.0,
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
                "bundle_axis_vs_rear_angle_deg": "",
                "hook_len_rel_err_max": 1.5,
            }
        ],
    )

    feature_schema = tmp_path / "features.yaml"
    feature_schema.write_text("feature_categories: {}\n", encoding="utf-8")
    campaign_config = {
        "kind": "generic_multi_run",
        "base_config": str(ROOT / "conf/sim_swim.yaml"),
        "base_overrides": {
            "time.duration_s": 0.1,
            "time.dt_star": 1.0e-4,
            "motor.torque_Nm": 2.5e-20,
            "motor.force_distribution": "root_torque_segment_couples",
        },
        "sweep": {
            "axes": {
                "attach_seed": {
                    "key": "seed.attach_seed",
                    "values": [0],
                    "ids": ["as000"],
                },
                "phase_seed": {
                    "key": "seed.phase_seed",
                    "values": [0],
                    "ids": ["ps000"],
                },
                "n_flagella": {
                    "key": "flagella.n_flagella",
                    "values": [1, 2],
                    "ids": ["nf01", "nf02"],
                },
            }
        },
        "dataset": {
            "dataset_id": "test_dataset",
            "feature_schema": str(feature_schema),
            "output_dir": str(dataset_dir),
            "timeseries_sampling": "all_steps",
        },
        "output": {
            "base_dir": str(run_dir),
            "timestamp_subdir": False,
        },
    }
    campaign_config_path = tmp_path / "campaign.yaml"
    campaign_config_path.write_text(yaml.safe_dump(campaign_config), encoding="utf-8")
    run_manifest = {
        "kind": "generic_multi_run",
        "output_root": str(run_dir),
        "conditions": [
            {
                "condition_id": "as000__ps000__nf01",
                "condition_label": "as=0, ps=0, nf=1",
                "output_dir": str(raw_a),
                "axis_values": {
                    "attach_seed": 0,
                    "phase_seed": 0,
                    "n_flagella": 1,
                },
                "config_overrides": {
                    "flagella": {"n_flagella": 1},
                    "seed": {"attach_seed": 0, "phase_seed": 0},
                },
            },
            {
                "condition_id": "as000__ps000__nf02",
                "condition_label": "as=0, ps=0, nf=2",
                "output_dir": str(raw_b),
                "axis_values": {
                    "attach_seed": 0,
                    "phase_seed": 0,
                    "n_flagella": 2,
                },
                "config_overrides": {
                    "flagella": {"n_flagella": 2},
                    "seed": {"attach_seed": 0, "phase_seed": 0},
                },
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
        campaign_config_path=campaign_config_path,
        run_manifest_path=run_manifest_path,
        run_dir=None,
        overwrite=True,
        cli_overrides=[
            "dataset.dataset_id=cli_dataset",
            f"dataset.output_dir={tmp_path / 'datasets/cli_dataset'}",
        ],
    )

    dataset_dir = tmp_path / "datasets/cli_dataset"
    assert out_dir == dataset_dir
    assert (dataset_dir / "dataset_manifest.json").is_file()
    dataset_manifest = json.loads(
        (dataset_dir / "dataset_manifest.json").read_text(encoding="utf-8")
    )
    assert dataset_manifest["dataset_id"] == "cli_dataset"
    assert dataset_manifest["cli_overrides"] == [
        "dataset.dataset_id=cli_dataset",
        f"dataset.output_dir={tmp_path / 'datasets/cli_dataset'}",
    ]
    assert dataset_manifest["effective_campaign_config"]["dataset"][
        "output_dir"
    ] == str(dataset_dir)
    assert (dataset_dir / "campaign_config_used.yaml").is_file()
    assert (dataset_dir / "feature_schema_used.yaml").read_text(
        encoding="utf-8"
    ) == "feature_categories: {}\n"
    assert (dataset_dir / "configs/nf01_as000_ps000.yaml").is_file()
    with (dataset_dir / "summary.csv").open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 2
    assert rows[0]["quality_class"] == "strict_pass"
    assert float(rows[0]["cell_displacement"]) == 1.0
    assert float(rows[0]["cell_straightness"]) == 1.0
    assert math.isnan(float(rows[0]["flagella_axis_alignment"]))
    assert math.isnan(float(rows[0]["cell_flagella_axis_angle"]))
    assert int(rows[0]["missing_value_count"]) >= 8
    assert rows[1]["quality_class"] == "relaxed_pass"
    assert rows[1]["hook_wrapped"] == "True"
    assert (dataset_dir / "qc_summary.csv").is_file()
    with (dataset_dir / "timeseries/nf01_as000_ps000.csv").open(
        encoding="utf-8", newline=""
    ) as handle:
        ts_rows = list(csv.DictReader(handle))
    assert len(ts_rows) == 2
    assert ts_rows[0]["dataset_id"] == "cli_dataset"
    assert math.isnan(float(ts_rows[-1]["flag_helix_axis_alignment_order"]))
    assert math.isnan(float(ts_rows[-1]["bundle_axis_vs_body_axis_angle_deg"]))


def test_plot_distributions_outputs_analysis_and_qc_artifacts(tmp_path: Path) -> None:
    dataset_dir = tmp_path / "datasets/test_dataset"
    dataset_dir.mkdir(parents=True)
    feature_schema = {
        "feature_categories": {
            "metadata": {"variables": ["sample_id", "n_flagella"]},
            "quality": {
                "variables": ["quality_class", "use_for_analysis"],
            },
            "cell_translation": {
                "ml_candidate": True,
                "variables": ["cell_mean_speed", "cell_straightness"],
            },
            "flagella_axis": {
                "ml_candidate": True,
                "variables": [
                    "flagella_axis_alignment",
                    "flagella_axis_pair_angle_mean",
                ],
            },
            "diagnostics": {
                "ml_candidate": False,
                "variables": ["hook_drift", "hook_wrapped", "first_fail_category"],
            },
        }
    }
    (dataset_dir / "feature_schema_used.yaml").write_text(
        yaml.safe_dump(feature_schema, sort_keys=False),
        encoding="utf-8",
    )
    (dataset_dir / "dataset_manifest.json").write_text(
        json.dumps({"dataset_id": "test_dataset"}),
        encoding="utf-8",
    )
    rows = [
        {
            "sample_id": "nf01_seed000",
            "dataset_id": "test_dataset",
            "n_flagella": 1,
            "quality_class": "strict_pass",
            "use_for_analysis": True,
            "shape_pass": True,
            "relaxed_pass": True,
            "review_required": False,
            "first_fail_category": "none",
            "missing_value_count": 2,
            "cell_mean_speed": 0.2,
            "cell_straightness": 0.9,
            "flagella_axis_alignment": "nan",
            "flagella_axis_pair_angle_mean": "nan",
            "hook_drift": 0.1,
            "hook_wrapped": False,
        },
        {
            "sample_id": "nf02_seed000",
            "dataset_id": "test_dataset",
            "n_flagella": 2,
            "quality_class": "relaxed_pass",
            "use_for_analysis": True,
            "shape_pass": False,
            "relaxed_pass": True,
            "review_required": True,
            "first_fail_category": "hook",
            "missing_value_count": 0,
            "cell_mean_speed": 0.4,
            "cell_straightness": 0.8,
            "flagella_axis_alignment": 0.95,
            "flagella_axis_pair_angle_mean": 10.0,
            "hook_drift": 1.0,
            "hook_wrapped": True,
        },
        {
            "sample_id": "nf03_seed000",
            "dataset_id": "test_dataset",
            "n_flagella": 3,
            "quality_class": "fail",
            "use_for_analysis": False,
            "shape_pass": False,
            "relaxed_pass": False,
            "review_required": True,
            "first_fail_category": "hook",
            "missing_value_count": 0,
            "cell_mean_speed": 0.1,
            "cell_straightness": 0.5,
            "flagella_axis_alignment": 0.7,
            "flagella_axis_pair_angle_mean": 30.0,
            "hook_drift": 2.0,
            "hook_wrapped": True,
        },
        {
            "sample_id": "nf06_seed000",
            "dataset_id": "test_dataset",
            "n_flagella": 6,
            "quality_class": "relaxed_pass",
            "use_for_analysis": True,
            "shape_pass": False,
            "relaxed_pass": True,
            "review_required": True,
            "first_fail_category": "hook",
            "missing_value_count": 0,
            "cell_mean_speed": 0.5,
            "cell_straightness": 0.7,
            "flagella_axis_alignment": 0.99,
            "flagella_axis_pair_angle_mean": 8.0,
            "hook_drift": 1.5,
            "hook_wrapped": True,
        },
    ]
    with (dataset_dir / "summary.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    result = plot_distributions.analyze_dataset(
        dataset_dir=dataset_dir,
        overwrite=False,
    )

    assert len(result["analysis_csvs"]) == 4
    assert all(path.is_file() for path in result["analysis_csvs"])
    assert result["distribution_plots"]
    assert result["qc_plots"]
    assert (
        dataset_dir / "plots/distributions/cell_translation_all_samples.png"
    ).is_file()
    assert (
        dataset_dir / "plots/distributions/cell_translation_use_for_analysis_true.png"
    ).is_file()
    assert (dataset_dir / "plots/qc/quality_class_by_n_flagella.png").is_file()

    nan_rows = list(
        csv.DictReader(
            (dataset_dir / "analysis/nan_summary.csv").open(
                encoding="utf-8",
                newline="",
            )
        )
    )
    target_nan = [
        row
        for row in nan_rows
        if row["feature"] == "flagella_axis_pair_angle_mean"
        and row["n_flagella"] != "all"
        and float(row["n_flagella"]) == 1.0
    ]
    assert target_nan
    assert int(target_nan[0]["nan_count"]) == 1

    quality_rows = list(
        csv.DictReader(
            (dataset_dir / "analysis/quality_summary.csv").open(
                encoding="utf-8",
                newline="",
            )
        )
    )
    assert any(
        row["metric"] == "use_for_analysis"
        and float(row["n_flagella"]) == 3.0
        and row["category"] == "False"
        and row["count"] == "1"
        for row in quality_rows
    )

    with pytest.raises(FileExistsError, match="--overwrite"):
        plot_distributions.analyze_dataset(
            dataset_dir=dataset_dir,
            overwrite=False,
        )

    plot_distributions.analyze_dataset(dataset_dir=dataset_dir, overwrite=True)
