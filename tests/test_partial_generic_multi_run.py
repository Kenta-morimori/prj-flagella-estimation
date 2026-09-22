from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
import yaml

from sim_swim.analysis.partial_generic_multi_run import export_completed_campaign
from sim_swim.analysis.flagella_count_behavior import save_state_archive
from sim_swim.analysis.phase2_replay import _load_inputs
from sim_swim.sim.core import SimulationState


def _run_summary(*, completed: bool) -> dict[str, object]:
    return {
        "execution": {
            "status": "completed" if completed else "running",
            "observed_final_t_s": 0.0 if completed else 0.00004,
            "expected_final_step_summary_t_s": 0.0,
        },
        "gates": {"shape_nonbody": {}, "shape_body": {}},
        "all_step_metrics": {},
    }


def test_export_completed_campaign_excludes_incomplete_conditions(
    tmp_path: Path,
) -> None:
    run_dir = tmp_path / "source"
    campaign_path = tmp_path / "campaign.yaml"
    campaign_path.write_text(
        yaml.safe_dump(
            {
                "base_config": "conf/sim_swim_2010.yaml",
                "base_overrides": {"time": {"duration_s": 0.0}},
                "sweep": {
                    "axes": {
                        "n_flagella": {
                            "key": "flagella.n_flagella",
                            "values": [1, 2],
                            "ids": ["nf01", "nf02"],
                        }
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    complete = run_dir / "nf01"
    complete.mkdir(parents=True)
    (complete / "run_summary.json").write_text(
        json.dumps(_run_summary(completed=True)), encoding="utf-8"
    )
    (complete / "state_archive.npz").write_bytes(b"archive")
    (complete / "trajectory.csv").write_text("t\n0\n", encoding="utf-8")

    incomplete = run_dir / "nf02"
    incomplete.mkdir()
    (incomplete / "run_summary.json").write_text(
        json.dumps(_run_summary(completed=False)), encoding="utf-8"
    )

    output_dir = export_completed_campaign(
        campaign_config=campaign_path,
        run_dir=run_dir,
        output_dir=tmp_path / "partial",
        overwrite=False,
    )

    manifest = json.loads((output_dir / "run_manifest.json").read_text())
    assert manifest["partial"] is True
    assert manifest["condition_order"] == ["nf01"]
    assert manifest["excluded_conditions"] == [
        {"condition_id": "nf02", "reason": "not_completed"}
    ]
    assert (output_dir / "summary.csv").is_file()
    assert (output_dir / "run.log").is_file()


def test_partial_checkpoint_replay_requires_two_explicit_opt_ins(
    tmp_path: Path,
) -> None:
    run_dir = tmp_path / "source"
    campaign_path = tmp_path / "campaign.yaml"
    campaign_path.write_text(
        yaml.safe_dump(
            {
                "base_config": "conf/sim_swim_2010.yaml",
                "base_overrides": {"time": {"duration_s": 0.0}},
                "sweep": {
                    "axes": {
                        "n_flagella": {
                            "key": "flagella.n_flagella",
                            "values": [1],
                            "ids": ["nf01"],
                        }
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    condition = run_dir / "nf01"
    condition.mkdir(parents=True)
    partial_summary = _run_summary(completed=False)
    partial_summary["execution"] = {"status": "partial"}
    (condition / "run_summary.json").write_text(
        json.dumps(partial_summary), encoding="utf-8"
    )
    save_state_archive(
        condition / "state_archive.partial.npz",
        [
            SimulationState(
                t=0.0,
                position_um=(0.0, 0.0, 0.0),
                quaternion=(0.0, 0.0, 0.0, 1.0),
                velocity_um_s=(0.0, 0.0, 0.0),
                omega_rad_s=(0.0, 0.0, 0.0),
                bead_positions_um=np.zeros((1, 3)),
            )
        ],
    )
    (condition / "trajectory.partial.csv").write_text("t\n0\n", encoding="utf-8")

    exported = export_completed_campaign(
        campaign_config=campaign_path,
        run_dir=run_dir,
        output_dir=tmp_path / "partial",
        overwrite=False,
        include_partial_checkpoint=True,
    )
    manifest = json.loads((exported / "run_manifest.json").read_text())
    assert manifest["partial_checkpoint_included"] is True
    assert manifest["conditions"][0]["partial_evidence"]["status"] == "partial"
    with pytest.raises(ValueError, match="--allow-partial"):
        _load_inputs(exported)
    rows, _, _ = _load_inputs(exported, allow_partial=True)
    assert rows[0]["condition_id"] == "nf01"
    archive_path = condition / "state_archive.partial.npz"
    missing_path = condition / "state_archive.partial.missing"
    archive_path.rename(missing_path)
    with pytest.raises(FileNotFoundError, match="Missing state archive"):
        _load_inputs(exported, allow_partial=True)
    missing_path.rename(archive_path)
    archive_path.write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="Invalid partial state archive"):
        _load_inputs(exported, allow_partial=True)
