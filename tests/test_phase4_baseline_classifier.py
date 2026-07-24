from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from flagella_estimation.phase4.baseline import (
    FEATURE_NAMES,
    classification_metrics,
    extract_clip_features,
    fit_nearest_centroid,
    predict_nearest_centroid,
)
from flagella_estimation.phase4.training import (
    Phase4BaselineConfig,
    load_baseline_config,
    train_baseline_classifier,
)


def _write_fixture_dataset(dataset_dir: Path, *, dataset_version: str = "v1") -> None:
    clips_dir = dataset_dir / "clips"
    clips_dir.mkdir(parents=True)
    split_rows: list[dict[str, object]] = []
    metadata_records = []
    split_names = ("train", "val", "test")
    for n_flagella in (1, 2, 3):
        for split_index, split in enumerate(split_names):
            clip_id = f"nf{n_flagella:02d}_{split}"
            group_key = f"phase2:{dataset_version}:{clip_id}"
            clip = np.zeros((13, 16, 16), dtype=np.uint8)
            for frame_index in range(13):
                offset = (frame_index + split_index) % 2
                for line_index in range(n_flagella):
                    row = 4 + line_index * 3
                    clip[frame_index, row, 3 + offset : 12 + offset] = (
                        80 + n_flagella * 50
                    )
            np.save(clips_dir / f"{clip_id}.npy", clip)
            split_rows.append(
                {
                    "clip_id": clip_id,
                    "sample_id": clip_id,
                    "group_key": group_key,
                    "split": split,
                    "n_flagella": n_flagella,
                }
            )
            metadata_records.append(
                {
                    "schema_version": "phase3_clip_metadata/v0",
                    "dataset_id": "phase4_baseline_fixture",
                    "source_video": {"source_kind": "phase2_pseudo"},
                    "processing_mode": "gt_passthrough",
                    "provenance": {
                        "dataset_version": dataset_version,
                        "model_id": "phase2_flagella_count_behavior_v1",
                        "render_id": "state_archive_numpy_v1",
                    },
                    "track": {"group_key": group_key},
                    "clip": {
                        "clip_id": clip_id,
                        "output_path": f"clips/{clip_id}.npy",
                        "frame_count": 13,
                        "duration_s": 0.5,
                        "window_policy": "non_overlap",
                    },
                    "normalization": {"crop_size_px": [16, 16]},
                    "labels": {
                        "n_flagella": n_flagella,
                        "label_source": "phase2_gt",
                    },
                    "frames": [{"clip_frame_index": index} for index in range(13)],
                    "qc": {"status": "pass", "exclusion_reason": None},
                }
            )

    (dataset_dir / "manifest.json").write_text(
        json.dumps(
            {
                "pipeline_name": "phase3_gt_passthrough",
                "schema_version": "phase3_clip_metadata/v0",
                "dataset_id": "phase4_baseline_fixture",
                "clip_count": len(metadata_records),
                "clip": {"duration_s": 0.5, "window_policy": "non_overlap"},
                "filters": {
                    "allowed_n_flagella": [1, 2, 3],
                    "require_use_for_ml_candidate": True,
                    "baseline_torque_Nm": 2.0e-20,
                },
            }
        ),
        encoding="utf-8",
    )
    (dataset_dir / "clip_metadata.jsonl").write_text(
        "".join(json.dumps(record) + "\n" for record in metadata_records),
        encoding="utf-8",
    )
    with (dataset_dir / "split_summary.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(split_rows[0]))
        writer.writeheader()
        writer.writerows(split_rows)


@pytest.mark.light
def test_baseline_features_and_metrics_are_deterministic() -> None:
    clip = np.zeros((3, 8, 8), dtype=np.uint8)
    clip[:, 2:4, 2:6] = 255
    features = extract_clip_features(clip)
    assert features.shape == (len(FEATURE_NAMES),)
    assert np.isfinite(features).all()

    train_features = np.asarray([[0.0, 0.0], [0.2, 0.1], [2.0, 2.0], [2.2, 2.1]])
    train_labels = np.asarray([1, 1, 2, 2])
    model = fit_nearest_centroid(train_features, train_labels)
    predicted = predict_nearest_centroid(model, train_features)
    metrics = classification_metrics(train_labels, predicted, model.classes)
    assert predicted.tolist() == train_labels.tolist()
    assert metrics == {
        "sample_count": 4,
        "accuracy": 1.0,
        "balanced_accuracy": 1.0,
        "macro_f1": 1.0,
    }


@pytest.mark.light
def test_train_baseline_classifier_writes_reproducible_artifacts(
    tmp_path: Path,
) -> None:
    dataset_dir = tmp_path / "phase3_dataset"
    output_dir = tmp_path / "phase4_baseline"
    _write_fixture_dataset(dataset_dir)

    result = train_baseline_classifier(
        Phase4BaselineConfig(dataset_dir=dataset_dir, output_dir=output_dir)
    )

    assert result == output_dir
    expected_files = {
        "confusion_matrix.csv",
        "manifest.json",
        "metrics.json",
        "model.npz",
        "predictions.csv",
        "run.log",
    }
    assert {path.name for path in output_dir.iterdir()} == expected_files
    metrics = json.loads((output_dir / "metrics.json").read_text(encoding="utf-8"))
    assert metrics["classes"] == [1, 2, 3]
    assert set(metrics["metrics_by_split"]) == {"train", "val", "test"}
    assert metrics["metrics_by_split"]["test"]["sample_count"] == 3

    prediction_rows = list(
        csv.DictReader(
            (output_dir / "predictions.csv").open("r", encoding="utf-8", newline="")
        )
    )
    assert len(prediction_rows) == 9
    assert {row["split"] for row in prediction_rows} == {"train", "val", "test"}
    assert len({row["group_key"] for row in prediction_rows}) == 9

    manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["freeze"] == {
        "allowed_n_flagella": [1, 2, 3],
        "clip_duration_s": 0.5,
        "dataset_version": "v1",
        "window_policy": "non_overlap",
    }
    assert manifest["dataset_audit"]["group_count"] == 9


@pytest.mark.light
def test_baseline_rejects_dataset_outside_freeze(tmp_path: Path) -> None:
    dataset_dir = tmp_path / "phase3_dataset"
    _write_fixture_dataset(dataset_dir, dataset_version="v2")

    with pytest.raises(ValueError, match="dataset_versions"):
        train_baseline_classifier(
            Phase4BaselineConfig(
                dataset_dir=dataset_dir,
                output_dir=tmp_path / "phase4_baseline",
            )
        )


@pytest.mark.light
def test_load_baseline_config_supports_key_value_overrides(tmp_path: Path) -> None:
    config_path = tmp_path / "baseline.yaml"
    config_path.write_text(
        "dataset_dir: input\noutput_dir: output\nseed: 1\n",
        encoding="utf-8",
    )
    cfg = load_baseline_config(
        config_path,
        [
            "dataset_dir=override_input",
            "output_dir=override_output",
            "freeze.allowed_n_flagella=[1, 2, 3]",
            "seed=7",
        ],
    )
    assert cfg.dataset_dir == Path("override_input")
    assert cfg.output_dir == Path("override_output")
    assert cfg.allowed_n_flagella == (1, 2, 3)
    assert cfg.seed == 7
