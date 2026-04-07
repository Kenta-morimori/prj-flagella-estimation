from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig


def _make_cfg(
    motor_torque_Nm: float = -1.0,
    hook_enabled: bool = True,
    n_flagella: int = 3,
    stub_mode: str = "full_flagella",
    duration_s: float = 5.0e-5,
    fd_eps_over_b: float = 1.0e-3,
) -> SimulationConfig:
    return SimulationConfig.from_dict(
        {
            "scale": {"b_um": 1.0, "bead_radius_a_over_b": 0.1},
            "body": {
                "prism": {
                    "n_prism": 3,
                    "dz_over_b": 0.5,
                    "radius_over_b": 0.5,
                    "axis": "x",
                },
                "length_total_um": 2.0,
            },
            "flagella": {
                "n_flagella": n_flagella,
                "placement_mode": "uniform",
                "init_mode": "paper_table1",
                "stub_mode": stub_mode,
                "n_beads_per_flagellum": 11,
                "discretization": {"ds_over_b": 0.58},
                "bond_L_over_b": 0.58,
                "length_over_b": 2.32,
                "helix_init": {"radius_over_b": 0.25, "pitch_over_b": 2.5},
            },
            "fluid": {"viscosity_Pa_s": 1.0e-3},
            "motor": {"torque_Nm": motor_torque_Nm, "reverse_n_flagella": 1},
            "potentials": {
                "spring": {"H_over_T_over_b": 10.0, "s": 0.1},
                "bend": {
                    "kb_over_T": 20.0,
                    "theta0_deg": {
                        "normal": 142.0,
                        "semicoiled": 90.0,
                        "curly1": 105.0,
                    },
                },
                "torsion": {
                    "kt_over_T": 10.0,
                    "fd_eps_over_b": fd_eps_over_b,
                    "phi0_deg": {
                        "normal": -60.0,
                        "semicoiled": 65.0,
                        "curly1": 120.0,
                    },
                },
                "spring_spring_repulsion": {
                    "A_ss_over_T": 1.0,
                    "a_ss_over_b": 0.2,
                    "cutoff_over_b": 0.2,
                },
            },
            "hook": {"enabled": hook_enabled, "threshold_deg": 90.0, "kb_over_T": 20.0},
            "run_tumble": {
                "run_tau": 20.0,
                "tumble_tau": 8.0,
                "semicoiled_tau": 4.0,
                "curly1_tau": 4.0,
            },
            "time": {"duration_s": duration_s, "dt_s": 1.0e-3},
            "output_sampling": {"out_all_steps_3d": True, "fps_out_2d": 25.0},
            "brownian": {
                "enabled": False,
                "temperature_K": 298.0,
                "method": "cholesky",
                "jitter": 1.0e-20,
            },
            "render": {
                "image_size_px": 128,
                "pixel_size_um": 0.203,
                "flagella_linewidth_px": 2.0,
                "render_flagella": True,
                "save_frames_3d": False,
                "follow_camera_3d": True,
                "view_range_um": 8.0,
                "timestamp_3d": False,
                "label_flagella": True,
                "follow_camera_2d": False,
                "save_frames_2d": False,
            },
            "seed": {"global_seed": 0},
            "output": {"base_dir": "outputs"},
        }
    )


def _run_and_load_step_summary(
    sim: Simulator, duration_s: float, summary_dir: Path
) -> list[dict[str, str]]:
    sim.run(duration_s, step_summary_dir=summary_dir)
    csv_path = summary_dir / "step_summary.csv"
    assert csv_path.is_file()
    with csv_path.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    assert rows
    return rows


def test_short_run_no_nan_inf() -> None:
    cfg = _make_cfg()
    sim = Simulator(cfg)
    states = sim.run(cfg.time.duration_s)

    assert len(states) >= 2

    arr_pos = np.array([s.position_um for s in states], dtype=float)
    arr_q = np.array([s.quaternion for s in states], dtype=float)
    arr_v = np.array([s.velocity_um_s for s in states], dtype=float)
    arr_w = np.array([s.omega_rad_s for s in states], dtype=float)

    assert np.isfinite(arr_pos).all()
    assert np.isfinite(arr_q).all()
    assert np.isfinite(arr_v).all()
    assert np.isfinite(arr_w).all()


def test_run_writes_step_summary_csv_without_projection_columns(tmp_path: Path) -> None:
    cfg = _make_cfg()
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, cfg.time.duration_s, tmp_path / "sim")

    first = rows[0]
    assert "step" in first
    assert "F_total_mean_all" in first
    assert "hook_count" in first
    assert "flag_intra_count" in first
    assert "flag_bend_err_mean_deg" in first
    assert "flag_torsion_err_mean_deg" in first
    assert "flag_state_changed" in first

    # projection は完全削除されていること
    assert "projection_body_rigid_enabled" not in first
    assert "projection_hook_length_enabled" not in first
    assert "projection_basal_link_direction_enabled" not in first
    assert "projection_flagella_chain_length_enabled" not in first
    assert "projection_flagella_template_enabled" not in first

    # 既存 diagnostics は維持
    for key in [
        "torsion_fd_eps_m",
        "torsion_fd_eps_over_b",
        "hook_len_mean_over_b",
        "flag_bond_len_mean_over_b",
        "flag_bond_rel_err_max",
        "local_attach_first_rel_err",
        "local_first_second_rel_err",
        "local_second_third_rel_err",
        "local_basal_bend_err_deg",
        "local_first_torsion_err_deg",
        "local_F_spring_attach_first",
        "local_F_motor_attach",
        "local_F_repulsion_basal_region",
    ]:
        assert key in first


def test_run_writes_initial_geometry_summary_json(tmp_path: Path) -> None:
    cfg = _make_cfg(motor_torque_Nm=0.0, hook_enabled=True, n_flagella=1)
    sim = Simulator(cfg)
    sim.run(cfg.time.duration_s, step_summary_dir=tmp_path / "sim")

    summary_path = tmp_path / "sim" / "initial_geometry_summary.json"
    assert summary_path.is_file()
    data = json.loads(summary_path.read_text(encoding="utf-8"))
    assert data["flagella"]["init_mode"] == "paper_table1"
    assert data["flagella"]["n_flagella"] == 1
    assert data["flagella"]["n_beads_per_flagellum"] == 11
    assert len(data["per_flagellum"]) == 1


def test_phase1_body_static_stability(tmp_path: Path) -> None:
    cfg = _make_cfg(
        motor_torque_Nm=0.0,
        hook_enabled=False,
        n_flagella=0,
        duration_s=0.01,
    )
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, cfg.time.duration_s, tmp_path / "phase1")

    assert len(rows) >= 10
    assert all(r["pos_all_finite"] in {"True", "true", "1"} for r in rows)
    assert all(int(r["hook_count"]) == 0 for r in rows)
    assert all(int(r["flag_intra_count"]) == 0 for r in rows)


def test_phase2_body_hook_static_stability_minimal_stub(tmp_path: Path) -> None:
    cfg = _make_cfg(
        motor_torque_Nm=0.0,
        hook_enabled=True,
        n_flagella=1,
        stub_mode="minimal_basal_stub",
        duration_s=0.01,
        fd_eps_over_b=1.0e-3,
    )
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, cfg.time.duration_s, tmp_path / "phase2")

    assert len(rows) >= 10
    assert all(r["pos_all_finite"] in {"True", "true", "1"} for r in rows)

    attach_rel = [float(r["local_attach_first_rel_err"]) for r in rows]
    first_second_rel = [float(r["local_first_second_rel_err"]) for r in rows]

    assert np.isfinite(np.array(attach_rel, dtype=float)).all()
    assert np.isfinite(np.array(first_second_rel, dtype=float)).all()

    # Phase2(静的)の暫定許容閾値: 5%
    assert max(attach_rel) < 0.05
    assert max(first_second_rel) < 0.05


def test_torsion_fd_eps_sweep_is_traceable_and_reduces_step0_torsion_force(
    tmp_path: Path,
) -> None:
    eps_values = [0.1, 0.001, 0.0001]
    torsion_force_step0: list[float] = []

    for eps in eps_values:
        cfg = _make_cfg(
            motor_torque_Nm=0.0,
            hook_enabled=True,
            n_flagella=1,
            stub_mode="full_flagella",
            duration_s=0.002,
            fd_eps_over_b=eps,
        )
        sim = Simulator(cfg)
        rows = _run_and_load_step_summary(
            sim, cfg.time.duration_s, tmp_path / f"sim_eps_{eps}"
        )
        first = rows[0]
        assert float(first["torsion_fd_eps_over_b"]) == pytest.approx(eps)
        torsion_force_step0.append(float(first["F_torsion_mean_flag"]))

    assert torsion_force_step0[1] < torsion_force_step0[0]
    assert torsion_force_step0[2] < torsion_force_step0[0]
