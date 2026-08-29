from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import pytest


SCRIPT = (
    Path(__file__).parents[1]
    / "scripts/03_dataset_building/stage_issue215_reference.py"
)
SPEC = importlib.util.spec_from_file_location("stage_issue215_reference", SCRIPT)
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def _campaign(root: Path, count: int = 36) -> None:
    conditions = []
    for index in range(count):
        condition_id = f"condition_{index:02d}"
        condition = root / "conditions" / condition_id
        condition.mkdir(parents=True)
        (condition / "run_summary.json").write_text(
            json.dumps(
                {
                    "execution": {"status": "completed"},
                    "gates": {
                        "finite": {
                            "status": "available",
                            "any_fail": False,
                            "final_pass": True,
                        },
                        "shape_nonbody": {
                            "status": "available",
                            "any_fail": False,
                            "final_pass": True,
                        },
                    },
                }
            ),
            encoding="utf-8",
        )
        conditions.append(
            {
                "condition_id": condition_id,
                "output_dir": str(condition),
                "axis_values": {"n_flagella": 1, "attach_seed": 0, "phase_seed": 0},
            }
        )
    (root / "run_manifest.json").write_text(
        json.dumps({"git": {"commit": "f2d0ce8"}, "conditions": conditions}),
        encoding="utf-8",
    )
    (root / "summary.csv").write_text("condition_id\n", encoding="utf-8")
    analysis = root / "analysis" / "motion_features"
    analysis.mkdir(parents=True)
    (analysis / "manifest.json").write_text("{}", encoding="utf-8")


@pytest.mark.light
def test_stage_issue215_reference_copies_only_compact_artifacts(tmp_path: Path) -> None:
    campaign, reference = tmp_path / "campaign", tmp_path / "reference"
    _campaign(campaign)
    MODULE.stage(campaign, reference)
    assert len(list(csv.DictReader((reference / "qc_summary.csv").open()))) == 36
    manifest = json.loads((reference / "reference_manifest.json").read_text())
    assert manifest["raw_archive_transferred"] is False
    assert manifest["run_summary_count"] == 36
    assert not list(reference.rglob("state_archive.npz"))


@pytest.mark.light
def test_stage_issue215_reference_rejects_incomplete_campaign(tmp_path: Path) -> None:
    campaign = tmp_path / "campaign"
    _campaign(campaign, count=35)
    with pytest.raises(ValueError, match="expected 36"):
        MODULE.stage(campaign, tmp_path / "reference")
