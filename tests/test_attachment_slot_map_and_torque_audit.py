from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

from sim_swim.analysis.attachment_slot_map import render_attachment_slot_map
from sim_swim.analysis.flagella_count_behavior import save_state_archive
from sim_swim.analysis.motor_torque_archive_audit import (
    _balance,
    _observed_match,
    _per_flag,
    audit_completed_archives,
)
from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    load_yaml,
    normalize_campaign_config,
)
from sim_swim.analysis.phase2_replay import _parse_args, _plot_cell
from sim_swim.analysis.torque_weight_replay import reconstructed_segment_weights
from sim_swim.render.render3d import _flagella_colors
from sim_swim.sim.core import SimulationState, Simulator
from sim_swim.sim.params import SimulationConfig

pytestmark = pytest.mark.light


def _state(positions_um: np.ndarray, t: float = 0.0) -> SimulationState:
    return SimulationState(
        t=t,
        position_um=(0.0, 0.0, 0.0),
        quaternion=(0.0, 0.0, 0.0, 1.0),
        velocity_um_s=(0.0, 0.0, 0.0),
        omega_rad_s=(0.0, 0.0, 0.0),
        bead_positions_um=positions_um,
    )


def test_replay_default_fixed_preserves_explicit_follow_and_does_not_recenter() -> None:
    assert _parse_args(["--input-dir", "/tmp/in"]).camera_3d == "fixed"
    assert (
        _parse_args(["--input-dir", "/tmp/in", "--camera-3d", "follow"]).camera_3d
        == "follow"
    )
    raw = json.loads(
        json.dumps(
            __import__("yaml").safe_load(
                Path("conf/sim_swim_2010_hex.yaml").read_text()
            )
        )
    )
    raw.setdefault("render", {}).pop("follow_camera_3d", None)
    cfg = SimulationConfig.from_dict(raw)
    assert cfg.render.follow_camera_3d is False
    simulator = Simulator(cfg)
    positions_um = simulator.model.positions_m * 1e6
    first = _state(positions_um)
    second = _state(positions_um + np.array([1.0, 0.0, 0.0]))
    second.position_um = (1.0, 0.0, 0.0)
    fig = plt.figure()
    axis = fig.add_subplot(projection="3d")
    _plot_cell(axis, st=first, cfg=cfg, rig=simulator.rig, title="first", fail_label="")
    before = axis.get_xlim()
    axis.cla()
    _plot_cell(
        axis, st=second, cfg=cfg, rig=simulator.rig, title="second", fail_label=""
    )
    assert axis.get_xlim() == before
    plt.close(fig)
    assert (
        cfg.with_overrides(
            {"render": {"follow_camera_3d": True}}
        ).render.follow_camera_3d
        is True
    )
    profile = SimulationConfig.from_dict(
        __import__("yaml").safe_load(Path("conf/sim_swim_2010.yaml").read_text())
    )
    assert profile.render.follow_camera_3d is True


def test_explicit_fixed_camera_is_static_and_validated() -> None:
    args = _parse_args(
        [
            "--input-dir",
            "/tmp/in",
            "--view",
            "3d",
            "--camera-3d",
            "fixed",
            "--view-range-mode",
            "explicit-fixed",
            "--camera-3d-center-um",
            "0",
            "0",
            "0",
            "--camera-3d-half-range-um",
            "2.8",
        ]
    )
    assert args.camera_3d_center_um == [0.0, 0.0, 0.0]
    assert args.camera_3d_half_range_um == 2.8
    with pytest.raises(SystemExit):
        _parse_args(
            [
                "--input-dir",
                "/tmp/in",
                "--view-range-mode",
                "explicit-fixed",
                "--camera-3d",
                "follow",
                "--camera-3d-center-um",
                "0",
                "0",
                "0",
                "--camera-3d-half-range-um",
                "2.8",
            ]
        )
    with pytest.raises(SystemExit):
        _parse_args(
            [
                "--input-dir",
                "/tmp/in",
                "--view-range-mode",
                "explicit-fixed",
                "--camera-3d-center-um",
                "0",
                "0",
                "0",
                "--camera-3d-half-range-um",
                "0",
            ]
        )

    cfg = SimulationConfig.from_dict(load_yaml(Path("conf/sim_swim_2010_hex.yaml")))
    simulator = Simulator(cfg)
    original = simulator.model.positions_m * 1e6
    fig = plt.figure()
    axis = fig.add_subplot(projection="3d")
    for offset in (0.0, 1.0):
        state = _state(original + np.array([offset, 0.0, 0.0]))
        state.position_um = (offset, 0.0, 0.0)
        axis.cla()
        _plot_cell(
            axis,
            st=state,
            cfg=cfg,
            rig=simulator.rig,
            title="fixed world origin",
            fail_label="",
            camera_center_um=np.zeros(3),
            view_range_um=2.8,
        )
        assert axis.get_xlim() == (-2.8, 2.8)
        assert axis.get_ylim() == (-2.8, 2.8)
        assert axis.get_zlim() == (-2.8, 2.8)
    plt.close(fig)


def test_slot_map_uses_archived_beads_and_3d_colors(tmp_path: Path) -> None:
    records = []
    for slots in ([0, 1, 2, 3], [0, 1, 2, 3, 4, 5]):
        n = len(slots)
        condition_id = f"nf{n:02d}__slots{''.join(map(str, slots))}"
        source = tmp_path / condition_id
        source.mkdir()
        ring = np.asarray(
            [[0.0, np.cos(i * np.pi / 3), np.sin(i * np.pi / 3)] for i in range(6)]
        )
        body = np.concatenate([ring + [float(layer - 2), 0, 0] for layer in range(5)])
        save_state_archive(source / "state_archive.npz", [_state(body)])
        records.append(
            {
                "condition_id": condition_id,
                "output_dir": str(source),
                "axis_values": {"attachment_slots": slots},
                "geometry": {
                    "actual": {
                        "body_slots_per_layer": 6,
                        "body_layers": 5,
                        "body_beads": 30,
                        "attachment_topology": [
                            {
                                "flag_id": i,
                                "layer": 2,
                                "slot": slot,
                                "body_bead_index": 12 + slot,
                            }
                            for i, slot in enumerate(slots)
                        ],
                    }
                },
            }
        )
    replay_input = tmp_path / "replay_input"
    replay_input.mkdir()
    (replay_input / "run_manifest.json").write_text(json.dumps({"conditions": records}))
    outputs = render_attachment_slot_map(replay_input, tmp_path / "map")
    assert outputs["attachment_slot_map"].is_file()
    manifest = json.loads(outputs["attachment_slot_map_manifest"].read_text())
    assert manifest["condition_count"] == 2
    assert manifest["view"] == "rear_to_front_along_positive_x"
    for n, row in zip((4, 6), manifest["conditions"], strict=True):
        assert [item["flag_id"] for item in row["mapping"]] == list(range(n))
        assert [item["slot"] for item in row["mapping"]] == list(range(n))
        assert [item["body_bead_index"] for item in row["mapping"]] == list(
            range(12, 12 + n)
        )
        assert [item["color_rgb"] for item in row["mapping"]] == [
            list(c) for c in _flagella_colors(n)
        ]
    assert manifest["conditions"][-1]["full_ring_rotation_equivalent"] is True


def test_tilted_flag_axial_only_leaves_transverse_torque_full_vector_cancels() -> None:
    body = (
        np.asarray(
            [[0, -1, -1], [0, 1, -1], [0, 1, 1], [0, -1, 1], [1, 0, 1]], dtype=float
        )
        * 1e-6
    )
    theta = np.linspace(0, 2 * np.pi, 11, endpoint=False)
    flag = (
        np.column_stack(
            [
                np.linspace(0, 2, 11),
                0.3 * np.cos(theta) + 0.12 * np.arange(11),
                0.25 * np.sin(theta),
            ]
        )
        * 1e-6
    )
    positions = np.vstack([body, flag])
    body_ids = np.arange(len(body))
    flag_ids = [np.arange(len(body), len(positions))]
    weights = [np.ones(10)]
    old, old_flags = _per_flag(positions, body_ids, flag_ids, weights, 2.5e-20, False)
    full, full_flags = _per_flag(positions, body_ids, flag_ids, weights, 2.5e-20, True)
    old_balance = _balance(positions, old, body_ids, flag_ids)
    full_balance = _balance(positions, full, body_ids, flag_ids)
    assert old_flags[0]["flag_transverse_torque_norm_Nm"] > 0
    assert old_balance["global_torque_residual_ratio"] > 0.01
    assert full_balance["global_torque_residual_ratio"] < 1e-10
    assert full_balance["global_force_residual_ratio"] < 1e-10
    assert np.linalg.norm(full_flags[0]["net_torque_Nm"]) < 1e-28


def test_recorded_mismatch_and_missing_archive_stop_counterfactual(
    tmp_path: Path,
) -> None:
    calculated = {
        "body_torque_Nm": [0.0, 0.0, 0.0],
        "flag_torque_Nm": [1e-20, 0.0, 0.0],
        "global_torque_residual_ratio": 0.1,
        "global_force_residual_ratio": 0.0,
    }
    observed = {
        **{
            f"diagnostic_motor_net_torque_{side}_{axis}_Nm": "0"
            for side in ("body", "flag")
            for axis in "xyz"
        },
        "diagnostic_motor_torque_balance_residual_ratio": "0.1",
        "diagnostic_motor_force_balance_residual_ratio": "0",
    }
    assert _observed_match(calculated, observed, 2.5e-20)["matched"] is False
    evaluation = tmp_path / "evaluation"
    replay = evaluation / "replay_input"
    replay.mkdir(parents=True)
    source = tmp_path / "source"
    source.mkdir()
    (replay / "run_manifest.json").write_text(
        json.dumps(
            {
                "base_config": "conf/sim_swim_2010_hex.yaml",
                "conditions": [{"condition_id": "missing", "output_dir": str(source)}],
            }
        )
    )
    with (evaluation / "summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["condition_id", "state_archive_sha256"]
        )
        writer.writeheader()
        writer.writerow(
            {
                "condition_id": "missing",
                "state_archive_sha256": hashlib.sha256(b"wrong").hexdigest(),
            }
        )
    with pytest.raises(FileNotFoundError, match="state_archive.npz"):
        audit_completed_archives(evaluation, tmp_path / "audit")
    (source / "state_archive.npz").write_bytes(b"not the recorded archive")
    (source / "run_summary.json").write_text('{"execution":{"status":"completed"}}')
    (source / "diagnostic_samples.csv").write_text("diagnostic_t_s\n0\n")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        audit_completed_archives(evaluation, tmp_path / "audit")


def test_completed_archive_match_unlocks_only_same_state_comparison(
    tmp_path: Path,
) -> None:
    campaign = normalize_campaign_config(
        load_yaml(
            Path(
                "conf/phase2_multi_run/2010_hex_project_long_duration_2s_issue245.yaml"
            )
        )
    )
    record = dict(build_campaign_conditions(campaign)[0])
    condition_id = record["condition_id"]
    source = tmp_path / "source"
    source.mkdir()
    record["output_dir"] = str(source)
    raw = load_yaml(Path("conf/sim_swim_2010_hex.yaml"))
    cfg = SimulationConfig.from_dict(raw).with_overrides(record["config_overrides"])
    simulator = Simulator(cfg)
    positions_um = simulator.model.positions_m * 1e6
    states = [_state(positions_um, 0.0), _state(positions_um, 0.01)]
    archive = source / "state_archive.npz"
    save_state_archive(archive, states)
    (source / "run_summary.json").write_text('{"execution":{"status":"completed"}}')
    weight = reconstructed_segment_weights(
        cfg.motor.torque_distribution_profile,
        len(simulator.model.flagella_indices[0]) - 1,
        times_s=np.asarray([0.01]),
        dt_s=cfg.dt_s,
        torque_Nm=cfg.motor.torque_Nm,
    )[0]
    forces, _ = _per_flag(
        simulator.model.positions_m,
        simulator.model.body_indices,
        simulator.model.flagella_indices,
        [weight],
        cfg.torque_for_forces_Nm,
        False,
    )
    balance = _balance(
        simulator.model.positions_m,
        forces,
        simulator.model.body_indices,
        simulator.model.flagella_indices,
    )
    observed = {
        "diagnostic_t_s": "0.01",
        "diagnostic_motor_torque_balance_residual_ratio": str(
            balance["global_torque_residual_ratio"]
        ),
        "diagnostic_motor_force_balance_residual_ratio": str(
            balance["global_force_residual_ratio"]
        ),
    }
    for side in ("body", "flag"):
        for axis, value in zip("xyz", balance[f"{side}_torque_Nm"], strict=True):
            observed[f"diagnostic_motor_net_torque_{side}_{axis}_Nm"] = str(value)
    with (source / "diagnostic_samples.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(observed))
        writer.writeheader()
        writer.writerow(observed)
    evaluation = tmp_path / "evaluation"
    replay = evaluation / "replay_input"
    replay.mkdir(parents=True)
    (replay / "run_manifest.json").write_text(
        json.dumps(
            {"base_config": "conf/sim_swim_2010_hex.yaml", "conditions": [record]}
        )
    )
    with (evaluation / "summary.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["condition_id", "state_archive_sha256"]
        )
        writer.writeheader()
        writer.writerow(
            {
                "condition_id": condition_id,
                "state_archive_sha256": hashlib.sha256(
                    archive.read_bytes()
                ).hexdigest(),
            }
        )
    result = audit_completed_archives(evaluation, tmp_path / "audit")
    states_audited = result["conditions"][0]["states"]
    assert states_audited["initial"]["full_vector_same_state"] is None
    assert states_audited["sampled_maximum"]["recorded_match"]["matched"] is True
    assert (
        states_audited["sampled_maximum"]["full_vector_same_state"][
            "global_torque_residual_ratio"
        ]
        < 1e-10
    )
