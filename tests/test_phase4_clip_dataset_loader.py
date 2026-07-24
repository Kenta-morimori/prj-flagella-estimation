from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from flagella_estimation.phase3.pipeline import Phase3Config, build_clip_dataset
from flagella_estimation.phase4.dataset import (
    audit_phase4_clip_dataset,
    load_phase3_common_clip_dataset,
)
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


def _write_phase2_summary(input_dataset: Path, rows: list[dict[str, str]]) -> None:
    input_dataset.mkdir(parents=True)
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
        writer.writerows(rows)


@pytest.mark.light
def test_phase4_loader_reads_phase3_common_clip_dataset(tmp_path: Path) -> None:
    input_dataset = tmp_path / "phase2_dataset"
    phase2_rows: list[dict[str, str]] = []
    for n_flagella in (1, 2, 3):
        for run_index in range(3):
            sample_id = f"nf{n_flagella:02d}_run{run_index}"
            raw_dir = tmp_path / "raw" / sample_id
            raw_dir.mkdir(parents=True)
            save_state_archive(
                raw_dir / "state_archive.npz", [_state(i) for i in range(26)]
            )
            phase2_rows.append(
                {
                    "sample_id": sample_id,
                    "n_flagella": str(n_flagella),
                    "torque_Nm": "2e-20",
                    "use_for_ml_candidate": "True",
                    "raw_dir": str(raw_dir),
                }
            )
    _write_phase2_summary(input_dataset, phase2_rows)

    dataset_dir = build_clip_dataset(
        Phase3Config(
            dataset_id="phase4_loader_fixture",
            input_dataset=input_dataset,
            output_dir=tmp_path / "phase3_clips",
            crop_size_px=32,
            max_per_class=3,
        )
    )

    samples = load_phase3_common_clip_dataset(dataset_dir)
    audit = audit_phase4_clip_dataset(samples)

    assert len(samples) == 18
    assert audit.sample_count == 18
    assert audit.group_count == 9
    assert audit.split_counts == {"train": 6, "val": 6, "test": 6}
    assert audit.class_counts == {1: 6, 2: 6, 3: 6}
    assert audit.split_class_counts == {
        ("train", 1): 2,
        ("train", 2): 2,
        ("train", 3): 2,
        ("val", 1): 2,
        ("val", 2): 2,
        ("val", 3): 2,
        ("test", 1): 2,
        ("test", 2): 2,
        ("test", 3): 2,
    }
    assert {sample.frame_count for sample in samples} == {13}
    assert {sample.frame_shape for sample in samples} == {(32, 32)}
    assert {sample.n_flagella for sample in samples} == {1, 2, 3}


@pytest.mark.light
def test_phase4_loader_rejects_group_key_leakage(tmp_path: Path) -> None:
    input_dataset = tmp_path / "phase2_dataset"
    raw_dir = tmp_path / "raw" / "nf01_run0"
    raw_dir.mkdir(parents=True)
    save_state_archive(raw_dir / "state_archive.npz", [_state(i) for i in range(26)])
    _write_phase2_summary(
        input_dataset,
        [
            {
                "sample_id": "nf01_run0",
                "n_flagella": "1",
                "torque_Nm": "2e-20",
                "use_for_ml_candidate": "True",
                "raw_dir": str(raw_dir),
            }
        ],
    )
    dataset_dir = build_clip_dataset(
        Phase3Config(
            dataset_id="phase4_loader_fixture",
            input_dataset=input_dataset,
            output_dir=tmp_path / "phase3_clips",
            crop_size_px=32,
        )
    )

    split_path = dataset_dir / "split_summary.csv"
    rows = list(csv.DictReader(split_path.open("r", encoding="utf-8", newline="")))
    rows[1]["split"] = "val"
    with split_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    samples = load_phase3_common_clip_dataset(dataset_dir)
    with pytest.raises(ValueError, match="group_key leakage"):
        audit_phase4_clip_dataset(samples)


@pytest.mark.light
def test_phase4_loader_rejects_split_rows_without_metadata(tmp_path: Path) -> None:
    input_dataset = tmp_path / "phase2_dataset"
    raw_dir = tmp_path / "raw" / "nf01_run0"
    raw_dir.mkdir(parents=True)
    save_state_archive(raw_dir / "state_archive.npz", [_state(i) for i in range(26)])
    _write_phase2_summary(
        input_dataset,
        [
            {
                "sample_id": "nf01_run0",
                "n_flagella": "1",
                "torque_Nm": "2e-20",
                "use_for_ml_candidate": "True",
                "raw_dir": str(raw_dir),
            }
        ],
    )
    dataset_dir = build_clip_dataset(
        Phase3Config(
            dataset_id="phase4_loader_fixture",
            input_dataset=input_dataset,
            output_dir=tmp_path / "phase3_clips",
            crop_size_px=32,
        )
    )

    split_path = dataset_dir / "split_summary.csv"
    rows = list(csv.DictReader(split_path.open("r", encoding="utf-8", newline="")))
    extra = dict(rows[0])
    extra["clip_id"] = "missing_metadata_clip"
    rows.append(extra)
    with split_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    with pytest.raises(ValueError, match="missing_from_metadata"):
        load_phase3_common_clip_dataset(dataset_dir)


@pytest.mark.light
def test_phase4_loader_rejects_clip_shape_mismatch(tmp_path: Path) -> None:
    input_dataset = tmp_path / "phase2_dataset"
    raw_dir = tmp_path / "raw" / "nf01_run0"
    raw_dir.mkdir(parents=True)
    save_state_archive(raw_dir / "state_archive.npz", [_state(i) for i in range(26)])
    _write_phase2_summary(
        input_dataset,
        [
            {
                "sample_id": "nf01_run0",
                "n_flagella": "1",
                "torque_Nm": "2e-20",
                "use_for_ml_candidate": "True",
                "raw_dir": str(raw_dir),
            }
        ],
    )
    dataset_dir = build_clip_dataset(
        Phase3Config(
            dataset_id="phase4_loader_fixture",
            input_dataset=input_dataset,
            output_dir=tmp_path / "phase3_clips",
            crop_size_px=32,
        )
    )

    first_record = json.loads(
        (dataset_dir / "clip_metadata.jsonl")
        .read_text(encoding="utf-8")
        .splitlines()[0]
    )
    np.save(
        Path(first_record["clip"]["output_path"]),
        np.zeros((12, 32, 32), dtype=np.uint8),
    )

    with pytest.raises(ValueError, match="frame_count mismatch"):
        load_phase3_common_clip_dataset(dataset_dir)
