from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from sim_swim.analysis.flagella_count_behavior import save_state_archive
from sim_swim.analysis.initial_geometry_plot import render_initial_geometry
from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    load_yaml,
    normalize_campaign_config,
)
from sim_swim.sim.core import SimulationState


ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = (
    ROOT
    / "conf/phase2_multi_run/2010_hex_project_body_flagella_contact_screen_issue257.yaml"
)


def _state(positions: np.ndarray, *, time: float = 0.0) -> SimulationState:
    return SimulationState(
        t=time,
        position_um=(0.0, 0.0, 0.0),
        quaternion=(0.0, 0.0, 0.0, 1.0),
        velocity_um_s=(0.0, 0.0, 0.0),
        omega_rad_s=(0.0, 0.0, 0.0),
        bead_positions_um=positions,
    )


def _input(tmp_path: Path, *, changed_off: bool = False) -> Path:
    campaign = normalize_campaign_config(load_yaml(CAMPAIGN))
    conditions = build_campaign_conditions(campaign)
    records = []
    for condition in conditions:
        values = condition["axis_values"]
        slots = values["attachment_slots"]
        theta = np.arange(6) * np.pi / 3
        ring = np.column_stack((np.cos(theta), np.sin(theta)))
        body = np.array([[layer * 0.25, y, z] for layer in range(5) for y, z in ring])
        flags = []
        topology = []
        for index, slot in enumerate(slots):
            point = body[12 + slot]
            flags.extend([point + [0.1 + step * 0.15, 0.0, 0.0] for step in range(3)])
            topology.append(
                {
                    "flag_id": index,
                    "body_bead_index": 12 + slot,
                    "layer": 2,
                    "slot": slot,
                }
            )
        positions = np.vstack((body, np.asarray(flags)))
        if changed_off and condition["condition_id"] == "nf02__slots01__bfoff":
            positions[-1, 0] += 0.01
        output = tmp_path / condition["condition_id"]
        archive = output / "state_archive.npz"
        save_state_archive(archive, [_state(positions)])
        records.append(
            {
                "condition_id": condition["condition_id"],
                "output_dir": str(output),
                "axis_values": values,
                "artifact_sha256": {
                    "state_archive.npz": hashlib.sha256(
                        archive.read_bytes()
                    ).hexdigest()
                },
                "geometry": {
                    "actual": {
                        "body_slots_per_layer": 6,
                        "body_layers": 5,
                        "body_beads": 30,
                        "flagellum_beads": [3] * len(slots),
                        "attachment_topology": topology,
                    }
                },
            }
        )
    replay_input = tmp_path / "replay_input"
    replay_input.mkdir()
    (replay_input / "run_manifest.json").write_text(
        json.dumps({"conditions": records}), encoding="utf-8"
    )
    return replay_input


def test_initial_geometry_all_patterns_and_identical_arm_dedup(tmp_path: Path) -> None:
    replay_input = _input(tmp_path)
    outputs = render_initial_geometry(
        replay_input, tmp_path / "initial_geometry", b_um=1.0
    )
    manifest = json.loads(outputs["initial_geometry_manifest"].read_text())
    assert manifest["condition_count"] == 26
    assert manifest["panel_count"] == 13
    assert len(manifest["sources"]) == 26
    assert {item["n_flagella"] for item in manifest["panels"]} == set(range(1, 7))
    assert len(outputs) == 8
    assert all(len(item["condition_ids"]) == 2 for item in manifest["panels"])
    for name, expected_hash in manifest["images"].items():
        image = outputs["initial_geometry_overview"].parent / name
        assert image.is_file()
        assert hashlib.sha256(image.read_bytes()).hexdigest() == expected_hash


def test_initial_geometry_separates_distinct_on_off_and_rejects_bad_archive(
    tmp_path: Path,
) -> None:
    replay_input = _input(tmp_path, changed_off=True)
    outputs = render_initial_geometry(replay_input, tmp_path / "figures", b_um=1.0)
    manifest = json.loads(outputs["initial_geometry_manifest"].read_text())
    assert manifest["panel_count"] == 14
    assert sum(item["attachment_slots"] == [0, 1] for item in manifest["panels"]) == 2
    records_path = replay_input / "run_manifest.json"
    records = json.loads(records_path.read_text())
    records["conditions"][0]["artifact_sha256"]["state_archive.npz"] = "bad"
    records_path.write_text(json.dumps(records), encoding="utf-8")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        render_initial_geometry(replay_input, tmp_path / "invalid", b_um=1.0)


def test_initial_geometry_rejects_missing_or_noninitial_archive(tmp_path: Path) -> None:
    replay_input = _input(tmp_path)
    records_path = replay_input / "run_manifest.json"
    manifest = json.loads(records_path.read_text())
    record = manifest["conditions"][0]
    archive = Path(record["output_dir"]) / "state_archive.npz"
    archive.unlink()
    with pytest.raises(FileNotFoundError):
        render_initial_geometry(replay_input, tmp_path / "missing", b_um=1.0)
    save_state_archive(archive, [_state(np.zeros((33, 3)), time=0.001)])
    record["artifact_sha256"]["state_archive.npz"] = hashlib.sha256(
        archive.read_bytes()
    ).hexdigest()
    records_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="t=0"):
        render_initial_geometry(replay_input, tmp_path / "late", b_um=1.0)
