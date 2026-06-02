from __future__ import annotations

import csv
from dataclasses import asdict
from pathlib import Path

from sim_swim.sim.core import Simulator
from sim_swim.sim.helix_retention_gate import (
    HELIX_RETENTION_BEND_ERR_MAX_DEG_LIMIT,
    HELIX_RETENTION_BOND_REL_ERR_MAX_LIMIT,
    HELIX_RETENTION_TORSION_ERR_MAX_DEG_LIMIT,
    summarize_single_flagellum_helix_retention,
)
from sim_swim.sim.params import SimulationConfig


def _make_cfg(
    *,
    motor_torque_Nm: float = 1.0e-18,
    duration_s: float = 1.0,
    dt_star: float | None = None,
    local_spring_scale: float | None = None,
    local_bend_scale: float | None = None,
    local_torsion_scale: float | None = None,
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
                "n_flagella": 1,
                "placement_mode": "uniform",
                "init_mode": "paper_table1",
                "stub_mode": "full_flagella",
                "n_beads_per_flagellum": 11,
                "discretization": {"ds_over_b": 0.58},
                "bond_L_over_b": 0.58,
                "length_over_b": 2.32,
                "helix_init": {"radius_over_b": 0.25, "pitch_over_b": 2.5},
            },
            "fluid": {"viscosity_Pa_s": 1.0e-3},
            "motor": {
                "torque_Nm": motor_torque_Nm,
                "reverse_n_flagella": 1,
                "enable_switching": False,
                "local_hook_scale": 8.0,
                "local_spring_scale": 5.0,
                "local_bend_scale": 4.0,
                "local_torsion_scale": 4.0,
            },
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
                    "fd_eps_over_b": 1.0e-3,
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
            "hook": {"enabled": True, "threshold_deg": 90.0, "kb_over_T": 20.0},
            "run_tumble": {
                "run_tau": 20.0,
                "tumble_tau": 8.0,
                "semicoiled_tau": 4.0,
                "curly1_tau": 4.0,
            },
            "time": {"duration_s": duration_s, "dt_s": 1.0e-3, "dt_star": dt_star},
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


def _with_motor_local_scales(
    cfg: SimulationConfig,
    *,
    local_hook_scale: float | None = None,
    local_spring_scale: float | None = None,
    local_bend_scale: float | None = None,
    local_torsion_scale: float | None = None,
) -> SimulationConfig:
    cfg_dict = asdict(cfg)
    if local_hook_scale is not None:
        cfg_dict["motor"]["local_hook_scale"] = local_hook_scale
    if local_spring_scale is not None:
        cfg_dict["motor"]["local_spring_scale"] = local_spring_scale
    if local_bend_scale is not None:
        cfg_dict["motor"]["local_bend_scale"] = local_bend_scale
    if local_torsion_scale is not None:
        cfg_dict["motor"]["local_torsion_scale"] = local_torsion_scale
    return SimulationConfig.from_dict(cfg_dict)


def _with_motor_force_distribution(
    cfg: SimulationConfig,
    force_distribution: str,
) -> SimulationConfig:
    cfg_dict = asdict(cfg)
    cfg_dict["motor"]["force_distribution"] = force_distribution
    return SimulationConfig.from_dict(cfg_dict)


def _run_step_summary(cfg: SimulationConfig, out_dir: Path) -> list[dict[str, str]]:
    sim = Simulator(cfg)
    sim.run(cfg.time.duration_s, step_summary_dir=out_dir)

    csv_path = out_dir / "step_summary.csv"
    with csv_path.open("r", encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def test_flag_state_is_fixed_when_switching_disabled(tmp_path: Path) -> None:
    cfg = _make_cfg()
    rows = _run_step_summary(cfg, tmp_path / "sim")

    assert len(rows) >= 1000
    assert all(row["flag_state_changed"] in {"False", "false", "0"} for row in rows)


def test_phase26_default_break_fails_helix_retention_gate(tmp_path: Path) -> None:
    """P2-6-005: P2-5 break representative remains a reproducible flag fail."""
    cfg = _make_cfg(motor_torque_Nm=4.0e-21, duration_s=0.05)
    rows = _run_step_summary(cfg, tmp_path / "phase26_default_break")

    summary = summarize_single_flagellum_helix_retention(rows, min_steps=10)

    assert summary["helix_retention_pass"] is False
    assert summary["first_fail_category"] == "flag"
    assert summary["first_fail_step"] == 7
    assert summary["max_flag_bond_rel_err"] > HELIX_RETENTION_BOND_REL_ERR_MAX_LIMIT
    assert summary["max_flag_bend_err_deg"] > HELIX_RETENTION_BEND_ERR_MAX_DEG_LIMIT
    assert (
        summary["max_flag_torsion_err_deg"] > HELIX_RETENTION_TORSION_ERR_MAX_DEG_LIMIT
    )


def test_phase26_small_dt_without_bend_scale_loses_rotation(tmp_path: Path) -> None:
    """P2-6-005: dt縮小だけでは形状維持と回転activityを両立しない。"""
    cfg = _make_cfg(motor_torque_Nm=4.0e-21, duration_s=0.05, dt_star=2.5e-4)
    rows = _run_step_summary(cfg, tmp_path / "phase26_small_dt_only")

    summary = summarize_single_flagellum_helix_retention(rows, min_steps=100)

    assert summary["helix_retention_pass"] is False
    assert summary["first_fail_category"] == "motor_no_rotation"
    assert summary["max_flag_bond_rel_err"] < HELIX_RETENTION_BOND_REL_ERR_MAX_LIMIT
    assert summary["max_flag_bend_err_deg"] < HELIX_RETENTION_BEND_ERR_MAX_DEG_LIMIT
    assert (
        summary["max_flag_torsion_err_deg"] < HELIX_RETENTION_TORSION_ERR_MAX_DEG_LIMIT
    )
    assert summary["median_abs_flag_phase_rate_hz"] < 1.0


def test_phase26_small_dt_and_bend_scale_retain_helix(tmp_path: Path) -> None:
    """P2-6-005: old representative keeps shape but does not spin the helix."""
    cfg = _make_cfg(motor_torque_Nm=4.0e-21, duration_s=0.05, dt_star=2.5e-4)
    cfg = _with_motor_local_scales(
        cfg,
        local_hook_scale=8.0,
        local_spring_scale=None,
        local_bend_scale=8.0,
        local_torsion_scale=None,
    )
    rows = _run_step_summary(cfg, tmp_path / "phase26_small_dt_bend8")

    summary = summarize_single_flagellum_helix_retention(rows, min_steps=100)

    assert summary["helix_retention_pass"] is False
    assert summary["first_fail_category"] == "motor_no_rotation"
    assert summary["median_abs_flag_phase_rate_hz"] > 1.0
    assert summary["median_abs_flag_helix_spin_rate_hz"] < 1.0
    assert summary["max_flag_bond_rel_err"] < HELIX_RETENTION_BOND_REL_ERR_MAX_LIMIT
    assert summary["max_flag_bend_err_deg"] < HELIX_RETENTION_BEND_ERR_MAX_DEG_LIMIT
    assert (
        summary["max_flag_torsion_err_deg"] < HELIX_RETENTION_TORSION_ERR_MAX_DEG_LIMIT
    )


def test_phase26_higher_torque_jitter_does_not_count_as_net_spin(
    tmp_path: Path,
) -> None:
    """P2-6-005: high instantaneous spin must not pass without net rotation."""
    cfg = _make_cfg(motor_torque_Nm=8.0e-21, duration_s=0.05, dt_star=1.25e-4)
    cfg = _with_motor_local_scales(
        cfg,
        local_hook_scale=8.0,
        local_spring_scale=5.0,
        local_bend_scale=8.0,
        local_torsion_scale=4.0,
    )
    rows = _run_step_summary(cfg, tmp_path / "phase26_spin_bend8")

    summary = summarize_single_flagellum_helix_retention(rows, min_steps=200)

    assert summary["helix_retention_pass"] is False
    assert summary["first_fail_category"] == "motor_no_rotation"
    assert summary["step_count"] >= 399
    assert summary["median_abs_flag_helix_spin_rate_hz"] > 1.0
    assert summary["net_abs_flag_helix_spin_revolutions"] < 0.01
    assert summary["flag_helix_spin_direction_consistency"] < 0.1
    assert summary["min_flag_helix_spin_fit_r2"] > 0.5
    assert summary["max_flag_bond_rel_err"] < HELIX_RETENTION_BOND_REL_ERR_MAX_LIMIT
    assert summary["max_flag_bend_err_deg"] < HELIX_RETENTION_BEND_ERR_MAX_DEG_LIMIT
    assert (
        summary["max_flag_torsion_err_deg"] < HELIX_RETENTION_TORSION_ERR_MAX_DEG_LIMIT
    )
    assert summary["max_hook_len_rel_err"] < 0.5


def test_phase26_distributed_torque_drives_net_helix_spin(tmp_path: Path) -> None:
    """P2-6-006: distributed torque drives monotonic net helix rotation."""
    cfg = _make_cfg(motor_torque_Nm=2.0e-20, duration_s=0.05, dt_star=1.0e-4)
    cfg = _with_motor_local_scales(
        cfg,
        local_hook_scale=1.0,
        local_spring_scale=1.2,
        local_bend_scale=1.0,
        local_torsion_scale=1.0,
    )
    cfg = _with_motor_force_distribution(cfg, "distributed_flagellum")
    rows = _run_step_summary(cfg, tmp_path / "phase26_distributed_spin")

    summary = summarize_single_flagellum_helix_retention(
        rows,
        min_steps=400,
        min_net_abs_spin_revolutions=0.1,
        min_direction_consistency=0.3,
    )

    assert summary["helix_retention_pass"] is True
    assert summary["first_fail_category"] == "none"
    assert summary["net_abs_flag_helix_spin_revolutions"] > 0.1
    assert summary["flag_helix_spin_direction_consistency"] > 0.9
    assert summary["max_hook_len_rel_err"] < 0.5
    assert summary["max_flag_bond_rel_err"] < HELIX_RETENTION_BOND_REL_ERR_MAX_LIMIT
    assert summary["max_flag_bend_err_deg"] < HELIX_RETENTION_BEND_ERR_MAX_DEG_LIMIT
    assert (
        summary["max_flag_torsion_err_deg"] < HELIX_RETENTION_TORSION_ERR_MAX_DEG_LIMIT
    )


def test_phase26_axial_torque_flux_probe_drives_net_helix_spin(
    tmp_path: Path,
) -> None:
    """P2-6-007: axial torque flux proxy transmits torque better than triplet."""
    cfg = _make_cfg(motor_torque_Nm=2.0e-20, duration_s=0.05, dt_star=1.0e-4)
    cfg = _with_motor_local_scales(
        cfg,
        local_hook_scale=1.0,
        local_spring_scale=1.2,
        local_bend_scale=1.0,
        local_torsion_scale=1.0,
    )
    cfg = _with_motor_force_distribution(cfg, "axial_torque_flux_probe")
    rows = _run_step_summary(cfg, tmp_path / "phase26_axial_torque_flux_probe")

    summary = summarize_single_flagellum_helix_retention(
        rows,
        min_steps=400,
        min_net_abs_spin_revolutions=0.1,
        min_direction_consistency=0.3,
    )

    assert summary["helix_retention_pass"] is True
    assert summary["first_fail_category"] == "none"
    assert summary["net_abs_flag_helix_spin_revolutions"] > 0.1
    assert summary["helix_to_root_net_rotation_ratio"] > 1.0
    assert summary["flag_helix_spin_direction_consistency"] > 0.9
    assert summary["max_hook_len_rel_err"] < 0.5
    assert summary["max_flag_bond_rel_err"] < HELIX_RETENTION_BOND_REL_ERR_MAX_LIMIT
    assert summary["max_flag_bend_err_deg"] < HELIX_RETENTION_BEND_ERR_MAX_DEG_LIMIT
    assert (
        summary["max_flag_torsion_err_deg"] < HELIX_RETENTION_TORSION_ERR_MAX_DEG_LIMIT
    )


def test_phase26_local_twist_transmission_probe_tracks_twist_state(
    tmp_path: Path,
) -> None:
    """P2-6-007: local twist state propagates root torque toward the tip."""
    cfg = _make_cfg(motor_torque_Nm=2.0e-20, duration_s=0.05, dt_star=1.0e-4)
    cfg = _with_motor_local_scales(
        cfg,
        local_hook_scale=1.0,
        local_spring_scale=1.2,
        local_bend_scale=1.0,
        local_torsion_scale=1.0,
    )
    cfg = _with_motor_force_distribution(cfg, "local_twist_transmission_probe")
    rows = _run_step_summary(cfg, tmp_path / "phase26_local_twist_probe")

    summary = summarize_single_flagellum_helix_retention(
        rows,
        min_steps=400,
        min_net_abs_spin_revolutions=0.1,
        min_direction_consistency=0.3,
    )
    last = rows[-1]

    assert summary["helix_retention_pass"] is True
    assert summary["first_fail_category"] == "none"
    assert summary["net_abs_flag_helix_spin_revolutions"] > 0.1
    assert summary["flag_helix_spin_direction_consistency"] > 0.9
    assert float(last["local_twist_root_orientation_deg"]) > 1.0
    assert float(last["local_twist_tip_orientation_deg"]) > 0.0
    assert float(last["local_twist_abs_max_deg"]) > 0.0
    assert float(last["local_twist_tip_activity_ratio"]) > 0.0
    assert summary["max_hook_len_rel_err"] < 0.5
    assert summary["max_flag_bond_rel_err"] < HELIX_RETENTION_BOND_REL_ERR_MAX_LIMIT
    assert summary["max_flag_bend_err_deg"] < HELIX_RETENTION_BEND_ERR_MAX_DEG_LIMIT
    assert (
        summary["max_flag_torsion_err_deg"] < HELIX_RETENTION_TORSION_ERR_MAX_DEG_LIMIT
    )
