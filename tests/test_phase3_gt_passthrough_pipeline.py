from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from flagella_estimation.phase3.metadata import build_gt_passthrough_metadata
from flagella_estimation.phase3.pipeline import (
    Phase3Config,
    build_clip_dataset,
    validate_training_candidate,
)
from flagella_estimation.phase3.replay import ReplayConfig, render_contact_sheet
from flagella_estimation.phase3.render import render_clip_array
from flagella_estimation.phase3.splits import (
    assign_grouped_splits,
    assert_no_group_leakage,
)
from flagella_estimation.phase3.windows import FrameWindow, generate_windows
from sim_swim.analysis.flagella_count_behavior import save_state_archive
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
        crop_size_px=32,
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
    assert "python" in manifest["environment"]
    assert manifest["sample_count"] == 1
    assert manifest["clip_count"] == 2
    assert manifest["dataset_version"] == "v1"
    assert manifest["dataset_revision"] is None
    assert manifest["dataset_summary"]["training_candidate_clip_count"] == 2
    assert len(metadata_lines) == 2
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
            crop_size_px=32,
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
            crop_size_px=32,
        )
    )

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
