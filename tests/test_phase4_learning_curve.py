from __future__ import annotations

import csv
from copy import deepcopy
import json
from pathlib import Path

import numpy as np
import pytest

from phase4_test_utils import write_phase4_fixture_dataset
from flagella_estimation.phase4.learning_curve import (
    GroupFeature,
    Phase4LearningCurveConfig,
    evaluate_grouped_learning_curve,
    load_learning_curve_config,
    validate_group_sets,
)


@pytest.mark.light
def test_grouped_learning_curve_writes_group_level_artifacts(
    tmp_path: Path,
) -> None:
    dataset_dir = tmp_path / "phase3_dataset"
    output_dir = tmp_path / "learning_curve"
    write_phase4_fixture_dataset(dataset_dir)

    result = evaluate_grouped_learning_curve(
        Phase4LearningCurveConfig(
            dataset_dir=dataset_dir,
            output_dir=output_dir,
            repeats=20,
            seed=11,
        )
    )

    assert result == output_dir
    assert {path.name for path in output_dir.iterdir()} == {
        "confusion.csv",
        "learning_curve.csv",
        "learning_curve_summary.csv",
        "manifest.json",
        "run.log",
    }
    curve_rows = list(
        csv.DictReader(
            (output_dir / "learning_curve.csv").open("r", encoding="utf-8", newline="")
        )
    )
    assert len(curve_rows) == 8
    assert {
        (row["train_groups_per_class"], row["train_group_count"]) for row in curve_rows
    } == {("1", "3")}
    assert (
        len(
            {
                row["selected_group_keys"]
                for row in curve_rows
                if row["train_groups_per_class"] == "1"
            }
        )
        == 8
    )
    for row in curve_rows:
        selected = set(json.loads(row["selected_group_keys"]))
        holdout = set(json.loads(row["holdout_group_keys"]))
        assert selected.isdisjoint(holdout)
    assert {row["evaluation_group_count"] for row in curve_rows} == {"3"}

    summary_rows = list(
        csv.DictReader(
            (output_dir / "learning_curve_summary.csv").open(
                "r", encoding="utf-8", newline=""
            )
        )
    )
    assert [
        (row["train_groups_per_class"], row["repeat_count"]) for row in summary_rows
    ] == [("1", "8")]
    assert "n_flagella_1_recall_mean" in summary_rows[0]
    assert "macro_f1_p02_5" in summary_rows[0]

    manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["learning_curve"]["unit"] == (
        "unique track.group_key per n_flagella class"
    )
    assert manifest["group_counts_by_split"] == {"test": 3, "train": 3, "val": 3}
    assert manifest["learning_curve"]["train_groups_per_class"] == [1]
    assert manifest["learning_curve"]["protected_split"] == "test"
    assert manifest["learning_curve"]["holdout_groups_per_class"] == 1


@pytest.mark.light
def test_grouped_learning_curve_excludes_diagnostic_clips(tmp_path: Path) -> None:
    dataset_dir = tmp_path / "phase3_dataset"
    output_dir = tmp_path / "learning_curve"
    write_phase4_fixture_dataset(dataset_dir)
    records = [
        json.loads(line)
        for line in (dataset_dir / "clip_metadata.jsonl")
        .read_text(encoding="utf-8")
        .splitlines()
    ]
    split_rows = list(
        csv.DictReader(
            (dataset_dir / "split_summary.csv").open("r", encoding="utf-8", newline="")
        )
    )
    for record, split_row in zip(records[:6], split_rows[:6], strict=True):
        diagnostic = deepcopy(record)
        clip_id = f"{record['clip']['clip_id']}_diagnostic"
        diagnostic["clip"]["clip_id"] = clip_id
        diagnostic["clip"]["output_path"] = f"clips/{clip_id}.npy"
        diagnostic["qc"]["training_candidate"] = False
        diagnostic["qc"]["diagnostic_only"] = True
        np.save(
            dataset_dir / diagnostic["clip"]["output_path"],
            np.load(dataset_dir / record["clip"]["output_path"]),
        )
        records.append(diagnostic)
        duplicate_row = dict(split_row)
        duplicate_row["clip_id"] = clip_id
        split_rows.append(duplicate_row)
    (dataset_dir / "clip_metadata.jsonl").write_text(
        "".join(json.dumps(record) + "\n" for record in records),
        encoding="utf-8",
    )
    with (dataset_dir / "split_summary.csv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(split_rows[0]))
        writer.writeheader()
        writer.writerows(split_rows)
    manifest = json.loads((dataset_dir / "manifest.json").read_text(encoding="utf-8"))
    manifest["clip_count"] = len(records)
    (dataset_dir / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")

    evaluate_grouped_learning_curve(
        Phase4LearningCurveConfig(
            dataset_dir=dataset_dir,
            output_dir=output_dir,
            repeats=4,
        )
    )

    curve_rows = list(
        csv.DictReader(
            (output_dir / "learning_curve.csv").open("r", encoding="utf-8", newline="")
        )
    )
    assert {row["train_clip_count"] for row in curve_rows} == {"3"}


@pytest.mark.light
def test_grouped_learning_curve_is_deterministic_for_seed(tmp_path: Path) -> None:
    dataset_dir = tmp_path / "phase3_dataset"
    write_phase4_fixture_dataset(dataset_dir)
    first_dir = tmp_path / "first"
    second_dir = tmp_path / "second"

    for output_dir in (first_dir, second_dir):
        evaluate_grouped_learning_curve(
            Phase4LearningCurveConfig(
                dataset_dir=dataset_dir,
                output_dir=output_dir,
                repeats=4,
                seed=7,
            )
        )

    assert (first_dir / "learning_curve.csv").read_text(encoding="utf-8") == (
        second_dir / "learning_curve.csv"
    ).read_text(encoding="utf-8")
    assert (first_dir / "learning_curve_summary.csv").read_text(encoding="utf-8") == (
        second_dir / "learning_curve_summary.csv"
    ).read_text(encoding="utf-8")


@pytest.mark.light
def test_grouped_learning_curve_rejects_dataset_outside_freeze(
    tmp_path: Path,
) -> None:
    dataset_dir = tmp_path / "phase3_dataset"
    write_phase4_fixture_dataset(dataset_dir, dataset_version="v2")

    with pytest.raises(ValueError, match="dataset_versions"):
        evaluate_grouped_learning_curve(
            Phase4LearningCurveConfig(
                dataset_dir=dataset_dir,
                output_dir=tmp_path / "learning_curve",
            )
        )


@pytest.mark.light
def test_group_set_validation_rejects_development_evaluation_leakage() -> None:
    group = GroupFeature(
        group_key="phase2:v1:run0",
        split="train",
        n_flagella=1,
        features=np.zeros(2),
        clip_count=1,
    )
    with pytest.raises(ValueError, match="group leakage"):
        validate_group_sets([group], [group], np.asarray([1]))


@pytest.mark.light
def test_load_learning_curve_config_supports_key_value_overrides(
    tmp_path: Path,
) -> None:
    config_path = tmp_path / "curve.yaml"
    config_path.write_text(
        "dataset_dir: input\noutput_dir: output\nseed: 1\n",
        encoding="utf-8",
    )
    cfg = load_learning_curve_config(
        config_path,
        [
            "dataset_dir=override_input",
            "learning_curve.train_groups_per_class=[1, 2]",
            "learning_curve.repeats=5",
            "freeze.warmup_s=0.5",
            "seed=13",
        ],
    )
    assert cfg.dataset_dir == Path("override_input")
    assert cfg.train_groups_per_class == (1, 2)
    assert cfg.repeats == 5
    assert cfg.warmup_s == pytest.approx(0.5)
    assert cfg.seed == 13
