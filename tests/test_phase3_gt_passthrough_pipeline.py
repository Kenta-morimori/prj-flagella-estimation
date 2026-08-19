from __future__ import annotations

import csv
from dataclasses import replace
import json
from pathlib import Path

import numpy as np
import pytest

from flagella_estimation.phase3.metadata import build_gt_passthrough_metadata
from flagella_estimation.phase3.feature_comparison import pixel_features
from flagella_estimation.phase3.pipeline import (
    Phase3Config,
    build_clip_dataset,
    load_config,
    select_samples,
    validate_training_candidate,
)
from flagella_estimation.phase3.replay import (
    ReplayConfig,
    _load_replay_clip,
    render_3d_2d_grid_mp4,
    render_contact_sheet,
)
from flagella_estimation.phase3.render import render_clip_array
from flagella_estimation.phase3.splits import (
    assign_grouped_splits,
    assert_no_group_leakage,
)
from flagella_estimation.phase3.windows import FrameWindow, generate_windows
from flagella_estimation.phase3.feature_comparison import (
    PIXEL_FEATURES,
    grouped_nearest_centroid_pixel_baseline,
)
from sim_swim.analysis.flagella_count_behavior import save_state_archive
from sim_swim.render.body2d import render_body_capsule_frame
from sim_swim.sim.core import SimulationState


def _state(index: int) -> SimulationState:
    t_s = index / 25.0
    x_um = index * 0.01
    beads = np.asarray(
        [
            [x_um - 0.2, -0.05, 0.0],
            [x_um, 0.0, 0.0],
            [x_um + 0.2, 0.05, 0.0],
        ],
        dtype=float,
    )
    return SimulationState(
        t=t_s,
        position_um=(x_um, 0.0, 0.0),
        quaternion=(0.0, 0.0, 0.0, 1.0),
        velocity_um_s=(1.0, 0.0, 0.0),
        omega_rad_s=(0.0, 0.0, 0.0),
        bead_positions_um=beads,
        flag_states=(),
        reverse_flagella=(),
    )


def _state_with_flagella_far_from_body(index: int) -> SimulationState:
    state = _state(index)
    beads = np.asarray(
        [
            [state.position_um[0] - 0.2, -0.05, 0.0],
            [state.position_um[0] + 0.2, 0.05, 0.0],
            [state.position_um[0] + 3.0, 2.0, 0.0],
            [state.position_um[0] + 3.4, 2.4, 0.0],
        ],
        dtype=float,
    )
    return SimulationState(
        t=state.t,
        position_um=state.position_um,
        quaternion=state.quaternion,
        velocity_um_s=state.velocity_um_s,
        omega_rad_s=state.omega_rad_s,
        bead_positions_um=beads,
        flag_states=(),
        reverse_flagella=(),
    )


def _raw_state(index: int) -> SimulationState:
    t_s = index * 0.0001
    return SimulationState(
        t=t_s,
        position_um=(t_s, 0.0, 0.0),
        quaternion=(0.0, 0.0, 0.0, 1.0),
        velocity_um_s=(1.0, 0.0, 0.0),
        omega_rad_s=(0.0, 0.0, 0.0),
        bead_positions_um=np.asarray(
            [
                [t_s - 0.2, -0.05, 0.0],
                [t_s, 0.0, 0.0],
                [t_s + 0.2, 0.05, 0.0],
            ],
            dtype=float,
        ),
        flag_states=(),
        reverse_flagella=(),
    )


@pytest.mark.light
def test_phase3_window_generation_defaults_to_0p5s_non_overlap() -> None:
    windows = generate_windows(
        source_frame_count=26,
        frame_rate_hz=25.0,
        duration_s=0.5,
        policy="non_overlap",
    )

    assert windows == [FrameWindow(0, 13), FrameWindow(13, 26)]


@pytest.mark.light
def test_phase3_window_generation_supports_0p25s_1p0s_and_overlap() -> None:
    short = generate_windows(
        source_frame_count=26,
        frame_rate_hz=25.0,
        duration_s=0.25,
        policy="non_overlap",
    )
    long = generate_windows(
        source_frame_count=26,
        frame_rate_hz=25.0,
        duration_s=1.0,
        policy="non_overlap",
    )
    overlap = generate_windows(
        source_frame_count=26,
        frame_rate_hz=25.0,
        duration_s=0.5,
        policy="overlap",
        overlap_stride_fraction=0.5,
    )

    assert [window.frame_count for window in short] == [7, 7, 7]
    assert long == [FrameWindow(0, 25)]
    assert overlap[:3] == [FrameWindow(0, 13), FrameWindow(6, 19), FrameWindow(12, 25)]


@pytest.mark.light
def test_grouped_nearest_centroid_pixel_baseline_holds_out_source_groups() -> None:
    rows = []
    for n_flagella in (1, 2, 3):
        for seed in (0, 1):
            value = float(n_flagella * 10 + seed * 0.1)
            rows.append(
                {
                    "clip_id": f"n{n_flagella}_s{seed}",
                    "group_key": f"source-n{n_flagella}-s{seed}",
                    "n_flagella": n_flagella,
                    "training_candidate": True,
                    **{feature: value for feature in PIXEL_FEATURES},
                }
            )

    predictions, summary = grouped_nearest_centroid_pixel_baseline(rows)

    assert len(predictions) == 6
    assert summary["group_holdout"] is True
    assert summary["eligible_group_count"] == 6
    assert summary["accuracy"] == pytest.approx(1.0)


@pytest.mark.light
def test_phase3_grouped_split_rejects_cross_split_group_key() -> None:
    with pytest.raises(ValueError, match="group_key leakage"):
        assert_no_group_leakage(
            [
                {"group_key": "phase2:v1:run-a", "split": "train"},
                {"group_key": "phase2:v1:run-a", "split": "val"},
            ]
        )


@pytest.mark.light
def test_phase3_grouped_split_can_stratify_by_n_flagella() -> None:
    group_labels = {
        f"phase2:v1:nf{n_flagella:02d}_run{run_index}": n_flagella
        for n_flagella in (1, 2, 3)
        for run_index in range(3)
    }

    assignments = assign_grouped_splits(
        group_labels.keys(),
        group_labels=group_labels,
    )

    for n_flagella in (1, 2, 3):
        label_splits = {
            assignments[group_key]
            for group_key, label in group_labels.items()
            if label == n_flagella
        }
        assert label_splits == {"train", "val", "test"}


@pytest.mark.light
def test_phase3_config_uses_output_sampling_and_render_keys(tmp_path: Path) -> None:
    config_path = tmp_path / "phase3.yaml"
    config_path.write_text(
        "\n".join(
            [
                "dataset_id: phase3_fixture",
                f"input_dataset: {tmp_path / 'input'}",
                f"output_dir: {tmp_path / 'out'}",
                "clip:",
                "  duration_s: 0.5",
                "  window_policy: non_overlap",
                "  overlap_stride_fraction: 0.5",
                "output_sampling:",
                "  fps_out: 20.0",
                "render:",
                "  image_size_px: 40",
                "  pixel_size_um: 0.2",
                "  body_length_um: 2.5",
                "  body_width_um: 0.8",
            ]
        ),
        encoding="utf-8",
    )

    cfg = load_config(config_path)

    assert cfg.fps_out == pytest.approx(20.0)
    assert cfg.image_size_px == 40
    assert cfg.pixel_size_um == pytest.approx(0.2)
    assert cfg.body_length_um == pytest.approx(2.5)
    assert cfg.body_width_um == pytest.approx(0.8)


@pytest.mark.light
def test_phase3_keeps_configured_diagnostic_class_out_of_ml_scope() -> None:
    cfg = Phase3Config(
        dataset_id="fixture",
        input_dataset=Path("input"),
        output_dir=Path("out"),
        allowed_n_flagella=(1, 2, 3),
        diagnostic_n_flagella=(4,),
        baseline_torque_Nm=2.5e-20,
    )
    rows = [
        {
            "sample_id": "nf03",
            "n_flagella": "3",
            "torque_Nm": "2.5e-20",
            "use_for_ml_candidate": "true",
        },
        {
            "sample_id": "nf04",
            "n_flagella": "4",
            "torque_Nm": "2.5e-20",
            "use_for_ml_candidate": "true",
        },
    ]

    assert [row["sample_id"] for row in select_samples(rows, cfg)] == ["nf03", "nf04"]
    assert validate_training_candidate(rows[1], cfg) == (
        False,
        "n_flagella_not_in_mvp_scope",
    )


@pytest.mark.light
def test_phase3_pixel_features_measure_rendered_silhouette() -> None:
    frames = np.full((2, 16, 16), 255, dtype=np.uint8)
    frames[:, 6:10, 3:13] = 60
    features = pixel_features(frames, 25.0)

    assert features["2d_area_px_mean"] == pytest.approx(40.0)
    assert features["2d_aspect_ratio_mean"] > 1.0
    assert features["2d_axis_angular_velocity_rms_rad_s"] == pytest.approx(0.0)


@pytest.mark.light
def test_phase3_gt_passthrough_metadata_matches_required_schema_fields(
    tmp_path: Path,
) -> None:
    states = [_state(i) for i in range(13)]
    clip_array, geometries = render_clip_array(
        states,
        image_size_px=32,
        pixel_size_um=0.1,
    )
    clip_path = tmp_path / "clip.npy"
    np.save(clip_path, clip_array)
    source_path = tmp_path / "state_archive.npz"
    source_path.write_bytes(b"fixture")

    metadata = build_gt_passthrough_metadata(
        dataset_id="phase3_fixture",
        source_video_id="nf01_as000_ps000",
        source_path=source_path,
        frame_rate_hz=25.0,
        source_frame_count=13,
        source_duration_s=13 / 25.0,
        run_id="nf01_as000_ps000",
        raw_run_dir=tmp_path,
        n_flagella=1,
        track_id="nf01_as000_ps000:gt_track_0000",
        group_key="phase2:v1:nf01_as000_ps000",
        clip_id="nf01_as000_ps000_c0000",
        clip_index=0,
        window=FrameWindow(0, 13),
        window_policy="non_overlap",
        output_path=clip_path,
        crop_size_px=32,
        pixel_size_um=0.1,
        frame_geometries=geometries,
    )

    assert metadata["schema_version"] == "phase3_clip_metadata/v0"
    assert metadata["processing_mode"] == "gt_passthrough"
    assert metadata["labels"] == {"n_flagella": 1, "label_source": "phase2_gt"}
    assert metadata["provenance"]["run_id"] == "nf01_as000_ps000"
    assert metadata["provenance"]["render_id"] == "body_capsule_orthographic_v1"
    assert metadata["provenance"]["dataset_revision"] is None
    assert metadata["track"]["group_key"] == "phase2:v1:nf01_as000_ps000"
    assert metadata["track"]["source_frame_end"] == 12
    assert metadata["track"]["t_end_s"] == 12 / 25.0
    assert metadata["clip"]["frame_count"] == 13
    assert metadata["clip"]["source_frame_end"] == 12
    assert metadata["clip"]["t_end_s"] == 12 / 25.0
    assert metadata["clip"]["time_band"] == "0.00-0.52s"
    assert metadata["qc"]["training_candidate"] is True
    assert metadata["qc"]["diagnostic_only"] is False
    assert len(metadata["frames"]) == 13


@pytest.mark.light
def test_phase3_render_uses_body_only_rigid_capsule() -> None:
    states = [_state_with_flagella_far_from_body(i) for i in range(3)]

    clip_array, geometries = render_clip_array(
        states,
        image_size_px=64,
        pixel_size_um=0.1,
        body_length_um=2.0,
        body_width_um=1.0,
    )

    assert clip_array.shape == (3, 64, 64)
    assert clip_array.dtype == np.uint8
    assert np.min(clip_array) < 255
    # The far-away flagella-like beads are ignored, so the distant upper-right
    # crop area stays background.
    assert np.all(clip_array[:, 48:, 48:] == 255)
    assert {geometry.body_length_px for geometry in geometries} == {20.0}
    assert {geometry.body_width_px for geometry in geometries} == {10.0}


@pytest.mark.light
def test_phase3_render_orthographic_capsule_foreshortening() -> None:
    base = _state_with_flagella_far_from_body(0)
    parallel = replace(base, quaternion=(0.0, 0.0, 0.0, 1.0))
    tilt_60_deg = replace(base, quaternion=(0.0, 0.5, 0.0, 3**0.5 / 2.0))
    end_on = replace(
        base,
        quaternion=(0.0, 2**0.5 / 2.0, 0.0, 2**0.5 / 2.0),
    )

    _, geometries = render_clip_array(
        [parallel, tilt_60_deg, end_on],
        image_size_px=64,
        pixel_size_um=0.1,
        body_length_um=2.0,
        body_width_um=1.0,
    )

    assert geometries[0].body_length_px == pytest.approx(20.0)
    assert geometries[1].body_length_px == pytest.approx(15.0)
    assert geometries[2].body_length_px == pytest.approx(10.0)
    assert geometries[0].bbox_xywh_px[2:] == pytest.approx((20.0, 10.0))
    assert geometries[2].bbox_xywh_px[2:] == pytest.approx((10.0, 10.0))
    assert geometries[2].body_axis_angle_rad is None


@pytest.mark.light
def test_phase3_training_freeze_rejects_torque_variation_for_mvp(
    tmp_path: Path,
) -> None:
    cfg = Phase3Config(
        dataset_id="phase3_fixture",
        input_dataset=tmp_path,
        output_dir=tmp_path / "out",
    )
    baseline = {
        "n_flagella": "2",
        "use_for_ml_candidate": "True",
        "torque_Nm": "2e-20",
    }
    varied = dict(baseline)
    varied["torque_Nm"] = "4e-20"
    nf4 = dict(baseline)
    nf4["n_flagella"] = "4"

    assert validate_training_candidate(baseline, cfg) == (True, None)
    assert validate_training_candidate(varied, cfg) == (
        False,
        "torque_variation_diagnostic_only",
    )
    assert validate_training_candidate(nf4, cfg) == (
        False,
        "n_flagella_not_in_mvp_scope",
    )
    cfg_with_diagnostics = Phase3Config(
        dataset_id="phase3_fixture",
        input_dataset=tmp_path,
        output_dir=tmp_path / "out_with_diagnostics",
        source_require_use_for_ml_candidate=False,
    )
    diagnostic_source_row = dict(baseline)
    diagnostic_source_row["use_for_ml_candidate"] = "False"
    assert validate_training_candidate(diagnostic_source_row, cfg_with_diagnostics) == (
        True,
        None,
    )


@pytest.mark.light
def test_phase3_pipeline_writes_clips_manifest_and_summaries(tmp_path: Path) -> None:
    input_dataset = tmp_path / "dataset"
    raw_dir = tmp_path / "raw" / "nf01_as000_ps000"
    raw_dir.mkdir(parents=True)
    save_state_archive(raw_dir / "state_archive.npz", [_state(i) for i in range(26)])
    input_dataset.mkdir()
    with (input_dataset / "summary.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sample_id",
                "n_flagella",
                "torque_Nm",
                "use_for_ml_candidate",
                "raw_dir",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "sample_id": "nf01_as000_ps000",
                "n_flagella": "1",
                "torque_Nm": "2e-20",
                "use_for_ml_candidate": "True",
                "raw_dir": str(raw_dir),
            }
        )

    output_dir = tmp_path / "out"
    cfg = Phase3Config(
        dataset_id="phase3_fixture",
        input_dataset=input_dataset,
        output_dir=output_dir,
        config_path=Path("conf/phase3/fixture.yaml"),
        cli_overrides=("clip.duration_s=0.5",),
        image_size_px=32,
    )
    result_dir = build_clip_dataset(cfg)

    assert result_dir == output_dir
    manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
    metadata_lines = (
        (output_dir / "clip_metadata.jsonl").read_text(encoding="utf-8").splitlines()
    )
    clip = np.load(output_dir / "clips" / "nf01_as000_ps000_c0000.npy")

    assert manifest["schema_version"] == "phase3_clip_metadata/v0"
    assert manifest["invocation"] == {
        "config_path": "conf/phase3/fixture.yaml",
        "cli_overrides": ["clip.duration_s=0.5"],
    }
    assert manifest["filters"]["max_per_class"] is None
    assert manifest["clip"] == {
        "duration_s": 0.5,
        "window_policy": "non_overlap",
        "overlap_stride_fraction": 0.5,
    }
    assert manifest["output_sampling"] == {"fps_out": 25.0}
    assert manifest["render"] == {
        "render_mode": "body_capsule_orthographic_v1",
        "render_id": manifest["render"]["render_id"],
        "projection": "orthographic",
        "body_deformation_rendered": False,
        "rendered_objects": ["body"],
        "excluded_objects": ["flagella"],
        "body_shape": "capsule",
        "body_length_definition": "end_to_end",
        "image_size_px": 32,
        "pixel_size_um": 0.1,
        "body_length_um": 2.0,
        "body_width_um": 1.0,
        "body_intensity": 60,
        "background_intensity": 255,
        "tracking_center": True,
    }
    assert "python" in manifest["environment"]
    assert manifest["sample_count"] == 1
    assert manifest["clip_count"] == 2
    assert manifest["dataset_version"] == "v1"
    assert manifest["dataset_revision"] is None
    assert manifest["dataset_summary"]["training_candidate_clip_count"] == 2
    assert len(metadata_lines) == 2
    render_id = manifest["render"]["render_id"]
    assert render_id.startswith("body_capsule_orthographic_v1-variant-")
    assert {json.loads(line)["provenance"]["render_id"] for line in metadata_lines} == {
        render_id
    }
    assert clip.shape == (13, 32, 32)
    assert clip.dtype == np.uint8
    assert (output_dir / "run.log").is_file()
    assert (output_dir / "split_summary.csv").is_file()
    assert (output_dir / "qc_summary.csv").is_file()


@pytest.mark.light
def test_phase3_pipeline_marks_first_fail_windows_diagnostic(tmp_path: Path) -> None:
    input_dataset = tmp_path / "dataset"
    raw_dir = tmp_path / "raw" / "nf03_as001_ps001"
    raw_dir.mkdir(parents=True)
    save_state_archive(raw_dir / "state_archive.npz", [_state(i) for i in range(76)])
    input_dataset.mkdir()
    with (input_dataset / "summary.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sample_id",
                "n_flagella",
                "torque_Nm",
                "use_for_ml_candidate",
                "shape_pass",
                "first_fail_t_s",
                "raw_dir",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "sample_id": "nf03_as001_ps001",
                "n_flagella": "3",
                "torque_Nm": "2e-20",
                "use_for_ml_candidate": "False",
                "shape_pass": "False",
                "first_fail_t_s": "1.78",
                "raw_dir": str(raw_dir),
            }
        )

    output_dir = tmp_path / "out"
    build_clip_dataset(
        Phase3Config(
            dataset_id="phase3_fixture",
            input_dataset=input_dataset,
            output_dir=output_dir,
            image_size_px=32,
            dataset_revision="r1",
            source_require_use_for_ml_candidate=False,
        )
    )
    records = [
        json.loads(line)
        for line in (output_dir / "clip_metadata.jsonl").read_text().splitlines()
    ]

    assert len(records) == 5
    assert [record["qc"]["training_candidate"] for record in records] == [
        True,
        True,
        True,
        False,
        False,
    ]
    assert records[3]["qc"]["exclusion_reason"] == "post_or_contains_first_fail"
    assert records[0]["provenance"]["dataset_revision"] == "r1"
    summary = list(csv.DictReader((output_dir / "dataset_summary.csv").open()))[0]
    assert summary["training_candidate_clip_count"] == "3"
    assert summary["diagnostic_clip_count"] == "2"


@pytest.mark.light
def test_phase3_replay_writes_contact_sheet_and_manifest(tmp_path: Path) -> None:
    input_dataset = tmp_path / "dataset"
    raw_dir = tmp_path / "raw" / "nf01_as000_ps000"
    raw_dir.mkdir(parents=True)
    save_state_archive(raw_dir / "state_archive.npz", [_state(i) for i in range(26)])
    input_dataset.mkdir()
    with (input_dataset / "summary.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sample_id",
                "n_flagella",
                "torque_Nm",
                "use_for_ml_candidate",
                "raw_dir",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "sample_id": "nf01_as000_ps000",
                "n_flagella": "1",
                "torque_Nm": "2e-20",
                "use_for_ml_candidate": "True",
                "raw_dir": str(raw_dir),
            }
        )
    dataset_dir = build_clip_dataset(
        Phase3Config(
            dataset_id="phase3_fixture",
            input_dataset=input_dataset,
            output_dir=tmp_path / "clips",
            image_size_px=32,
            body_length_um=2.5,
            body_width_um=0.8,
            body_intensity=90,
        )
    )

    first_record = json.loads(
        (dataset_dir / "clip_metadata.jsonl")
        .read_text(encoding="utf-8")
        .splitlines()[0]
    )
    loaded_clip = _load_replay_clip(
        first_record,
        ReplayConfig(dataset_dir=dataset_dir),
    )
    assert loaded_clip["body_render_cfg"].body_length_um == pytest.approx(2.5)
    assert loaded_clip["body_render_cfg"].body_width_um == pytest.approx(0.8)
    assert loaded_clip["body_render_cfg"].body_intensity == 90
    replay_frame, _ = render_body_capsule_frame(
        loaded_clip["states"][0], loaded_clip["body_render_cfg"]
    )
    saved_clip = np.load(dataset_dir / first_record["clip"]["output_path"])
    assert np.array_equal(replay_frame, saved_clip[0])

    replay_dir = render_contact_sheet(
        ReplayConfig(
            dataset_dir=dataset_dir,
            output_dir=tmp_path / "replay",
            n_flagella=1,
            training_candidate=True,
            max_clips=1,
            frames_per_clip=2,
        )
    )

    assert (replay_dir / "contact_sheet.png").is_file()
    manifest = json.loads((replay_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["clip_count"] == 1
    assert manifest["filters"]["n_flagella"] == 1
    assert "commit" in manifest["git"]
    assert {"python", "platform", "numpy", "opencv", "matplotlib"} <= set(
        manifest["environment"]
    )


@pytest.mark.light
def test_phase3_replay_writes_3d_2d_mp4_grid_and_manifest(tmp_path: Path) -> None:
    input_dataset = tmp_path / "dataset"
    raw_dir = tmp_path / "raw" / "nf01_as000_ps000"
    raw_dir.mkdir(parents=True)
    save_state_archive(raw_dir / "state_archive.npz", [_state(i) for i in range(26)])
    input_dataset.mkdir()
    with (input_dataset / "summary.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sample_id",
                "n_flagella",
                "torque_Nm",
                "use_for_ml_candidate",
                "raw_dir",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "sample_id": "nf01_as000_ps000",
                "n_flagella": "1",
                "torque_Nm": "2e-20",
                "use_for_ml_candidate": "True",
                "raw_dir": str(raw_dir),
            }
        )
    dataset_dir = build_clip_dataset(
        Phase3Config(
            dataset_id="phase3_fixture",
            input_dataset=input_dataset,
            output_dir=tmp_path / "clips",
            image_size_px=32,
        )
    )

    replay_dir = render_3d_2d_grid_mp4(
        ReplayConfig(
            dataset_dir=dataset_dir,
            output_dir=tmp_path / "replay",
            n_flagella=1,
            training_candidate=True,
            max_clips=1,
            frames_per_clip=2,
            panel_layout="3d+2d",
            panel_width_px=320,
            panel_height_px=160,
        )
    )

    assert (replay_dir / "3d_2d_grid.mp4").is_file()
    manifest = json.loads((replay_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["pipeline_name"] == "phase3_clip_replay_3d_2d_grid"
    assert manifest["clip_count"] == 1
    assert manifest["panel_layout"] == "3d+2d"
    assert manifest["video_count"] == 1
    assert manifest["videos"][0]["fps"] == 25.0
    assert manifest["videos"][0]["frame_count"] == 13
    assert manifest["video"]["fps"] == 25.0
    assert manifest["outputs"]["mp4_grid"].endswith("3d_2d_grid.mp4")
    assert manifest["outputs"]["mp4_grids"][0].endswith("3d_2d_grid.mp4")
    assert "commit" in manifest["git"]
    assert {"python", "platform", "numpy", "opencv", "matplotlib"} <= set(
        manifest["environment"]
    )


@pytest.mark.light
def test_phase3_replay_splits_all_matching_clips_across_mp4_grids(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    records = [{"clip": {"clip_id": f"clip-{index}"}} for index in range(5)]
    monkeypatch.setattr(
        "flagella_estimation.phase3.replay._select_records", lambda _cfg: records
    )
    monkeypatch.setattr(
        "flagella_estimation.phase3.replay._load_replay_clip",
        lambda record, *_args: {
            "record": record,
            "states": [object()],
            "body_render_cfg": object(),
            "sim_cfg": object(),
            "rig": object(),
        },
    )
    monkeypatch.setattr(
        "flagella_estimation.phase3.replay._render_clip_pair_panel",
        lambda *args, **kwargs: np.zeros((16, 16, 3), dtype=np.uint8),
    )
    monkeypatch.setattr(
        "flagella_estimation.phase3.replay._resolve_replay_fps", lambda _clips: 25.0
    )

    replay_dir = render_3d_2d_grid_mp4(
        ReplayConfig(
            dataset_dir=tmp_path / "dataset",
            output_dir=tmp_path / "replay",
            clips_per_video=2,
            panel_width_px=16,
            panel_height_px=16,
        )
    )

    manifest = json.loads((replay_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["clip_count"] == 5
    assert manifest["video_count"] == 3
    assert [video["clip_count"] for video in manifest["videos"]] == [2, 2, 1]
    assert [Path(path).name for path in manifest["outputs"]["mp4_grids"]] == [
        "3d_2d_grid_001.mp4",
        "3d_2d_grid_002.mp4",
        "3d_2d_grid_003.mp4",
    ]


@pytest.mark.light
def test_phase3_replay_maps_clip_indices_to_sampled_frames(tmp_path: Path) -> None:
    input_dataset = tmp_path / "dataset"
    raw_dir = tmp_path / "raw" / "nf01_as000_ps000"
    raw_dir.mkdir(parents=True)
    save_state_archive(
        raw_dir / "state_archive.npz", [_raw_state(i) for i in range(6001)]
    )
    input_dataset.mkdir()
    with (input_dataset / "summary.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "sample_id",
                "n_flagella",
                "torque_Nm",
                "use_for_ml_candidate",
                "raw_dir",
            ],
        )
        writer.writeheader()
        writer.writerow(
            {
                "sample_id": "nf01_as000_ps000",
                "n_flagella": "1",
                "torque_Nm": "2e-20",
                "use_for_ml_candidate": "True",
                "raw_dir": str(raw_dir),
            }
        )
    dataset_dir = build_clip_dataset(
        Phase3Config(
            dataset_id="phase3_fixture",
            input_dataset=input_dataset,
            output_dir=tmp_path / "clips",
            image_size_px=32,
        )
    )
    record = json.loads(
        (dataset_dir / "clip_metadata.jsonl")
        .read_text(encoding="utf-8")
        .splitlines()[0]
    )

    replay_clip = _load_replay_clip(
        record,
        ReplayConfig(dataset_dir=dataset_dir, max_clips=1),
    )

    assert [round(state.t, 2) for state in replay_clip["states"]] == [
        round(frame["t_s"], 2) for frame in record["frames"]
    ]
    assert replay_clip["states"][-1].t == pytest.approx(0.48)
