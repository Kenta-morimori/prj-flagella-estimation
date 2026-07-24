from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import yaml


def write_phase4_fixture_dataset(
    dataset_dir: Path, *, dataset_version: str = "v1"
) -> None:
    clips_dir = dataset_dir / "clips"
    clips_dir.mkdir(parents=True)
    split_rows: list[dict[str, object]] = []
    metadata_records = []
    source_records = []
    for n_flagella in (1, 2, 3):
        for split_index, split in enumerate(("train", "val", "test")):
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
                    "dataset_id": "phase4_learning_curve_fixture",
                    "source_video": {"source_kind": "phase2_pseudo"},
                    "processing_mode": "gt_passthrough",
                    "provenance": {
                        "dataset_version": dataset_version,
                        "model_id": "phase2_flagella_count_behavior_v1",
                        "render_id": "state_archive_numpy_v1",
                        "run_id": clip_id,
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
            source_records.append({"sample_id": clip_id, "n_flagella": n_flagella})

    source_dir = dataset_dir.parent / f"{dataset_dir.name}_phase2_source"
    _write_phase2_source_fixture(
        source_dir, source_records, dataset_version=dataset_version
    )

    (dataset_dir / "manifest.json").write_text(
        json.dumps(
            {
                "pipeline_name": "phase3_gt_passthrough",
                "schema_version": "phase3_clip_metadata/v0",
                "dataset_id": "phase4_learning_curve_fixture",
                "input_dataset": str(source_dir),
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


def _write_phase2_source_fixture(
    source_dir: Path,
    records: list[dict[str, object]],
    *,
    dataset_version: str,
) -> None:
    configs_dir = source_dir / "configs"
    configs_dir.mkdir(parents=True)
    samples = []
    for record in records:
        sample_id = str(record["sample_id"])
        n_flagella = int(record["n_flagella"])
        config_path = configs_dir / f"{sample_id}.yaml"
        config_path.write_text(
            yaml.safe_dump(
                {
                    "motor": {
                        "torque_Nm": 2.0e-20,
                        "enable_switching": False,
                    },
                    "brownian": {"enabled": False},
                    "flagella": {"n_flagella": n_flagella},
                },
                sort_keys=True,
            ),
            encoding="utf-8",
        )
        samples.append(
            {
                "sample_id": sample_id,
                "n_flagella": n_flagella,
                "status": "completed",
                "config_path": str(config_path),
            }
        )
    (source_dir / "dataset_manifest.json").write_text(
        json.dumps(
            {
                "dataset_id": dataset_version,
                "effective_campaign_config": {
                    "metadata": {
                        "dataset_version": dataset_version,
                        "model_id": "flag_spring2p25_body2p5_candidate",
                    }
                },
                "samples": samples,
            }
        ),
        encoding="utf-8",
    )
