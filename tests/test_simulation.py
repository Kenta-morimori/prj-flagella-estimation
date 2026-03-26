from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest

from sim_swim.dynamics.engine import StepDiagnostics
from sim_swim.sim.collapse_diagnostics import (
    COLLAPSE_DIAGNOSTICS_COLUMNS,
    COLLAPSE_SUMMARY_COLUMNS,
)
from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig


def _triplet_angle_rad(positions_m: np.ndarray, triplet: np.ndarray) -> float:
    i, j, k = (int(triplet[0]), int(triplet[1]), int(triplet[2]))
    u = positions_m[i] - positions_m[j]
    v = positions_m[k] - positions_m[j]
    nu = max(float(np.linalg.norm(u)), 1e-18)
    nv = max(float(np.linalg.norm(v)), 1e-18)
    c = float(np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0))
    return float(np.arccos(c))


def _body_spring_rows(sim: Simulator) -> np.ndarray:
    pairs = sim.model.spring_pairs
    bi = sim.model.bead_is_body[pairs[:, 0]]
    bj = sim.model.bead_is_body[pairs[:, 1]]
    return np.where(bi & bj)[0]


def _body_bending_rows(sim: Simulator) -> np.ndarray:
    return np.where(sim.model.bending_flag_ids < 0)[0]


def _measure_body_constraint_drift(sim: Simulator, n_steps: int) -> tuple[float, float]:
    model = sim.model
    body_spring_rows = _body_spring_rows(sim)
    body_bending_rows = _body_bending_rows(sim)
    theta0, _ = sim.engine._state_angles_rad()
    theta0_body = theta0[body_bending_rows].copy()
    spring_l0 = model.spring_rest_lengths_m[body_spring_rows].copy()
    spring_pairs = model.spring_pairs[body_spring_rows]

    max_stretch_ratio = 0.0
    max_angle_error_rad = 0.0
    for _ in range(max(1, int(n_steps))):
        sim.engine.step(sim.config.dt_star)
        pos = model.positions_m

        spring_l = np.linalg.norm(
            pos[spring_pairs[:, 1]] - pos[spring_pairs[:, 0]],
            axis=1,
        )
        stretch_ratio = np.abs(spring_l - spring_l0) / np.maximum(spring_l0, 1e-30)
        max_stretch_ratio = max(max_stretch_ratio, float(np.max(stretch_ratio)))

        theta = np.asarray(
            [
                _triplet_angle_rad(pos, model.bending_triplets[row])
                for row in body_bending_rows
            ],
            dtype=float,
        )
        angle_error_rad = np.abs(theta - theta0_body)
        max_angle_error_rad = max(max_angle_error_rad, float(np.max(angle_error_rad)))

    return max_stretch_ratio, max_angle_error_rad


def _make_cfg(
    motor_torque_Nm: float = -1.0,
    hook_enabled: bool = True,
    n_flagella: int = 3,
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
                    "phi0_deg": {"normal": -60.0, "semicoiled": 65.0, "curly1": 120.0},
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
            "time": {"duration_s": 5.0e-5, "dt_s": 1.0e-3},
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


def _assert_hook_stats_over_all_steps(
    rows: list[dict[str, str]], b_um: float, n_flagella: int
) -> None:
    target_um = 0.25 * b_um
    mean_tol_um = 0.05 * target_um
    min_limit_um = 0.85 * target_um
    max_limit_um = 1.15 * target_um

    for row in rows:
        assert int(row["hook_count"]) == n_flagella
        hook_mean = float(row["hook_len_mean_um"])
        hook_min = float(row["hook_len_min_um"])
        hook_max = float(row["hook_len_max_um"])
        assert abs(hook_mean - target_um) <= mean_tol_um
        assert hook_min >= min_limit_um
        assert hook_max <= max_limit_um


def _assert_flag_intra_stats_over_all_steps(
    rows: list[dict[str, str]], b_um: float
) -> None:
    target_um = 0.58 * b_um
    mean_tol_um = 0.02 * target_um
    min_limit_um = 0.95 * target_um
    max_limit_um = 1.05 * target_um

    for row in rows:
        assert int(row["flag_intra_count"]) > 0
        intra_mean = float(row["flag_intra_len_mean_um"])
        intra_min = float(row["flag_intra_len_min_um"])
        intra_max = float(row["flag_intra_len_max_um"])
        assert abs(intra_mean - target_um) <= mean_tol_um
        assert intra_min >= min_limit_um
        assert intra_max <= max_limit_um


def _assert_flag_angle_error_over_all_steps(
    rows: list[dict[str, str]],
    bend_mean_deg_limit: float,
    bend_max_deg_limit: float,
    torsion_mean_deg_limit: float,
    torsion_max_deg_limit: float,
) -> None:
    for row in rows:
        bend_mean = float(row["flag_bend_err_mean_deg"])
        bend_max = float(row["flag_bend_err_max_deg"])
        torsion_mean = float(row["flag_torsion_err_mean_deg"])
        torsion_max = float(row["flag_torsion_err_max_deg"])
        assert bend_mean <= bend_mean_deg_limit
        assert bend_max <= bend_max_deg_limit
        assert torsion_mean <= torsion_mean_deg_limit
        assert torsion_max <= torsion_max_deg_limit


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


def test_run_writes_step_summary_csv(tmp_path: Path) -> None:
    cfg = _make_cfg()
    sim = Simulator(cfg)
    sim.run(cfg.time.duration_s, step_summary_dir=tmp_path / "sim")

    step_csv = tmp_path / "sim" / "step_summary.csv"
    step_full_csv = tmp_path / "sim" / "step_summary_full.csv"
    assert step_csv.is_file()
    assert not step_full_csv.exists()

    with step_csv.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    assert rows
    first = rows[0]
    assert "step" in first
    assert "F_total_mean_all" in first
    assert "hook_count" in first
    assert "flag_intra_count" in first
    assert "flag_bend_err_mean_deg" in first
    assert "flag_torsion_err_mean_deg" in first
    assert "flag_state_changed" in first
    assert first["brownian_enabled"] in {"False", "false", "0"}


@pytest.mark.parametrize(
    ("motor_torque", "hook_enabled"),
    [
        (0.0, False),
        (0.0, True),
        (-1.0, True),
    ],
)
def test_body_flag_bond_metric_is_numeric_across_modes(
    tmp_path: Path, motor_torque: float, hook_enabled: bool
) -> None:
    cfg = _make_cfg(motor_torque_Nm=motor_torque, hook_enabled=hook_enabled)
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, cfg.time.duration_s, tmp_path / "sim")

    first = rows[0]
    assert first["pos_all_finite"] in {"True", "true", "1"}
    hook_len = float(first["hook_len_mean_um"])
    assert np.isfinite(hook_len)


@pytest.mark.parametrize("hook_enabled", [False, True])
def test_zero_torque_initial_steps_do_not_scatter(hook_enabled: bool) -> None:
    cfg = _make_cfg(motor_torque_Nm=0.0, hook_enabled=hook_enabled)
    sim = Simulator(cfg)

    mean_disp_um: list[float] = []
    for _ in range(5):
        diag = sim.engine.step(cfg.dt_star)
        disp = np.linalg.norm(diag.positions_after_m - diag.positions_before_m, axis=1)
        mean_disp_um.append(float(np.mean(disp) * 1e6))

    assert max(mean_disp_um) < 0.7


def test_hook_bond_stats_columns_and_expected_values(tmp_path: Path) -> None:
    cfg = _make_cfg(motor_torque_Nm=0.0, hook_enabled=True, n_flagella=3)
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, cfg.time.duration_s, tmp_path / "sim")

    first = rows[0]
    assert int(first["hook_count"]) == cfg.flagella.n_flagella
    mean_len = float(first["hook_len_mean_um"])
    min_len = float(first["hook_len_min_um"])
    max_len = float(first["hook_len_max_um"])
    target_um = 0.25 * cfg.scale.b_um
    tol_um = 0.05 * target_um
    assert abs(mean_len - target_um) <= tol_um
    assert abs(min_len - target_um) <= tol_um
    assert abs(max_len - target_um) <= tol_um


def test_body_constraints_match_initial_reference() -> None:
    cfg = _make_cfg(motor_torque_Nm=0.0, hook_enabled=False, n_flagella=0)
    sim = Simulator(cfg)
    model = sim.model

    body_spring_rows = _body_spring_rows(sim)
    spring_pairs = model.spring_pairs[body_spring_rows]
    spring_l = np.linalg.norm(
        model.positions_m[spring_pairs[:, 1]] - model.positions_m[spring_pairs[:, 0]],
        axis=1,
    )
    spring_l0 = model.spring_rest_lengths_m[body_spring_rows]
    rel_err = np.max(np.abs(spring_l - spring_l0) / np.maximum(spring_l0, 1e-30))
    assert rel_err < 1e-12

    body_bending_rows = _body_bending_rows(sim)
    theta = np.asarray(
        [
            _triplet_angle_rad(model.positions_m, model.bending_triplets[row])
            for row in body_bending_rows
        ],
        dtype=float,
    )
    theta0, _ = sim.engine._state_angles_rad()
    theta0_body = theta0[body_bending_rows]
    max_theta_err = float(np.max(np.abs(theta - theta0_body)))
    assert max_theta_err < 1e-12


def test_body_constraints_nearly_invariant_after_one_step_without_external_forces() -> (
    None
):
    cfg = _make_cfg(motor_torque_Nm=0.0, hook_enabled=False, n_flagella=0)
    sim = Simulator(cfg)
    model = sim.model

    body_spring_rows = _body_spring_rows(sim)
    body_bending_rows = _body_bending_rows(sim)
    theta0, _ = sim.engine._state_angles_rad()
    theta0_body = theta0[body_bending_rows].copy()
    spring_l0 = model.spring_rest_lengths_m[body_spring_rows].copy()

    sim.engine.step(cfg.dt_star)
    pos = model.positions_m

    spring_pairs = model.spring_pairs[body_spring_rows]
    spring_l = np.linalg.norm(pos[spring_pairs[:, 1]] - pos[spring_pairs[:, 0]], axis=1)
    rel_err = np.max(np.abs(spring_l - spring_l0) / np.maximum(spring_l0, 1e-30))
    assert rel_err < 1e-4

    theta = np.asarray(
        [
            _triplet_angle_rad(pos, model.bending_triplets[row])
            for row in body_bending_rows
        ],
        dtype=float,
    )
    max_theta_err = float(np.max(np.abs(theta - theta0_body)))
    assert max_theta_err < 1e-3


def test_body_rigidity_is_preserved_over_multiple_steps_without_external_forces() -> (
    None
):
    cfg = _make_cfg(motor_torque_Nm=0.0, hook_enabled=False, n_flagella=0)
    sim = Simulator(cfg)
    max_stretch_ratio, max_angle_error_rad = _measure_body_constraint_drift(
        sim, n_steps=1000
    )

    assert max_stretch_ratio < 1e-4, (
        f"body distance constraint drifted: max_stretch_ratio={max_stretch_ratio:.6e}"
    )
    assert max_angle_error_rad < 1e-3, (
        f"body angle constraint drifted: max_angle_error_rad={max_angle_error_rad:.6e}"
    )


def test_motor_on_1e19_stability_and_body_constraints() -> None:
    cfg = _make_cfg(motor_torque_Nm=1.0e-19, hook_enabled=True, n_flagella=1)
    sim = Simulator(cfg)
    max_stretch_ratio, max_angle_error_rad = _measure_body_constraint_drift(
        sim, n_steps=2000
    )

    assert np.isfinite(sim.model.positions_m).all()
    assert max_stretch_ratio < 5e-4, (
        f"motor-on body distance drifted: max_stretch_ratio={max_stretch_ratio:.6e}"
    )
    assert max_angle_error_rad < 2e-3, (
        f"motor-on body angle drifted: max_angle_error_rad={max_angle_error_rad:.6e}"
    )


def test_hook_length_multi_step_motor_off(tmp_path: Path) -> None:
    cfg = _make_cfg(
        motor_torque_Nm=0.0, hook_enabled=True, n_flagella=1
    ).with_overrides({"time": {"duration_s": 2.0}})
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, 2.0, tmp_path / "sim")

    assert len(rows) >= 2000
    _assert_hook_stats_over_all_steps(rows, cfg.scale.b_um, cfg.flagella.n_flagella)


def test_hook_length_multi_step_motor_on(tmp_path: Path) -> None:
    cfg = _make_cfg(
        motor_torque_Nm=1.0e-18,
        hook_enabled=True,
        n_flagella=1,
    ).with_overrides({"time": {"duration_s": 2.0}})
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, 2.0, tmp_path / "sim")

    assert len(rows) >= 2000
    _assert_hook_stats_over_all_steps(rows, cfg.scale.b_um, cfg.flagella.n_flagella)


def test_flag_intra_length_multi_step_motor_off(tmp_path: Path) -> None:
    cfg = _make_cfg(
        motor_torque_Nm=0.0, hook_enabled=True, n_flagella=1
    ).with_overrides({"time": {"duration_s": 2.0}})
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, 2.0, tmp_path / "sim")

    assert len(rows) >= 2000
    _assert_flag_intra_stats_over_all_steps(rows, cfg.scale.b_um)


def test_flag_intra_length_multi_step_motor_on(tmp_path: Path) -> None:
    cfg = _make_cfg(
        motor_torque_Nm=1.0e-18,
        hook_enabled=True,
        n_flagella=1,
    ).with_overrides({"time": {"duration_s": 2.0}})
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, 2.0, tmp_path / "sim")

    assert len(rows) >= 2000
    _assert_flag_intra_stats_over_all_steps(rows, cfg.scale.b_um)
    _assert_flag_angle_error_over_all_steps(
        rows,
        bend_mean_deg_limit=6.0,
        bend_max_deg_limit=20.0,
        torsion_mean_deg_limit=6.0,
        torsion_max_deg_limit=25.0,
    )


def test_template_projection_is_skipped_when_disabled(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = _make_cfg(
        motor_torque_Nm=0.0,
        hook_enabled=True,
        n_flagella=1,
    ).with_overrides(
        {
            "projection": {"enable_flagella_template_projection": False},
        }
    )
    sim = Simulator(cfg)

    calls = {"count": 0}

    def _count_calls(positions_m: np.ndarray) -> np.ndarray:
        calls["count"] += 1
        return positions_m

    monkeypatch.setattr(sim.engine, "_project_flagella_template", _count_calls)
    sim.engine._project_hook_and_flag_bonds(sim.model.positions_m.copy())
    assert calls["count"] == 0


def test_chain_length_projection_is_used_when_template_disabled(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = _make_cfg(
        motor_torque_Nm=0.0,
        hook_enabled=True,
        n_flagella=1,
    ).with_overrides(
        {
            "projection": {
                "enable_flagella_template_projection": False,
                "enable_flagella_chain_length_projection_when_template_off": True,
            },
        }
    )
    sim = Simulator(cfg)

    calls = {"count": 0}

    def _count_calls(positions_m: np.ndarray, iterations: int) -> np.ndarray:
        calls["count"] += 1
        return positions_m

    monkeypatch.setattr(sim.engine, "_project_flagella_chain_lengths", _count_calls)
    sim.engine._project_hook_and_flag_bonds(sim.model.positions_m.copy())
    assert calls["count"] == sim.engine.constraint_projection_iters


def test_template_projection_is_used_when_enabled(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = _make_cfg(
        motor_torque_Nm=0.0,
        hook_enabled=True,
        n_flagella=1,
    ).with_overrides(
        {
            "projection": {
                "enable_flagella_template_projection": True,
                "enable_flagella_chain_length_projection_when_template_off": True,
            },
        }
    )
    sim = Simulator(cfg)

    template_calls = {"count": 0}
    chain_calls = {"count": 0}

    def _count_template_calls(positions_m: np.ndarray) -> np.ndarray:
        template_calls["count"] += 1
        return positions_m

    def _count_chain_calls(positions_m: np.ndarray, iterations: int) -> np.ndarray:
        chain_calls["count"] += 1
        return positions_m

    monkeypatch.setattr(sim.engine, "_project_flagella_template", _count_template_calls)
    monkeypatch.setattr(
        sim.engine, "_project_flagella_chain_lengths", _count_chain_calls
    )
    sim.engine._project_hook_and_flag_bonds(sim.model.positions_m.copy())

    assert template_calls["count"] == sim.engine.constraint_projection_iters
    assert chain_calls["count"] == 0


def test_engine_repulsion_parameters_follow_config_overrides() -> None:
    cfg = _make_cfg(
        motor_torque_Nm=1.0e-18,
        hook_enabled=True,
        n_flagella=1,
    ).with_overrides(
        {
            "potentials": {
                "spring_spring_repulsion": {
                    "A_ss_over_T": 2.5,
                    "cutoff_over_b": 0.35,
                }
            }
        }
    )
    sim = Simulator(cfg)

    torque = max(cfg.torque_for_forces_Nm, 1e-30)
    assert sim.engine.repulsion_A == pytest.approx(2.5 * torque)
    assert sim.engine.repulsion_cutoff_m == pytest.approx(0.35 * cfg.b_m)


def test_engine_stiffness_scales_follow_config_overrides() -> None:
    cfg = _make_cfg(
        motor_torque_Nm=1.0e-18,
        hook_enabled=True,
        n_flagella=1,
    ).with_overrides(
        {
            "stiffness_scales": {
                "body": 220.0,
                "flag_bend": 330.0,
                "flag_torsion": 340.0,
            }
        }
    )
    sim = Simulator(cfg)

    assert sim.engine.body_stiffness_scale == pytest.approx(220.0)
    assert sim.engine.flag_bend_stiffness_scale == pytest.approx(330.0)
    assert sim.engine.flag_torsion_stiffness_scale == pytest.approx(340.0)


def _dummy_diag(positions_m: np.ndarray) -> StepDiagnostics:
    z = np.zeros_like(positions_m)
    return StepDiagnostics(
        dt_star=1.0e-3,
        dt_s=1.0e-3,
        positions_before_m=positions_m.copy(),
        positions_after_m=positions_m.copy(),
        spring_forces=z,
        bend_forces=z,
        torsion_forces=z,
        hook_forces=z,
        repulsion_forces=z,
        motor_forces=z,
        total_forces=z,
        brownian_enabled=False,
        brownian_disp_m=z,
        motor_degenerate_axis_count=0,
        motor_split_rank_deficient_count=0,
        motor_bond_length_clipped_count=0,
    )


def test_collapse_diagnostics_files_and_columns_are_written(tmp_path: Path) -> None:
    cfg = _make_cfg(
        motor_torque_Nm=0.0, hook_enabled=True, n_flagella=1
    ).with_overrides(
        {
            "time": {"duration_s": 0.005},
            "diagnostics": {"collapse": {"enabled": True, "write_every_step": True}},
        }
    )
    sim = Simulator(cfg)
    sim_dir = tmp_path / "sim"
    sim.run(0.005, step_summary_dir=sim_dir)

    diag_csv = tmp_path / "collapse_diagnostics.csv"
    summary_csv = tmp_path / "collapse_summary.csv"
    assert diag_csv.is_file()
    assert summary_csv.is_file()

    with diag_csv.open("r", encoding="utf-8", newline="") as f:
        diag_rows = list(csv.DictReader(f))
    assert diag_rows
    assert list(diag_rows[0].keys()) == COLLAPSE_DIAGNOSTICS_COLUMNS
    expected_steps = int(round(0.005 / cfg.dt_s))
    assert len(diag_rows) == expected_steps * cfg.flagella.n_flagella

    with summary_csv.open("r", encoding="utf-8", newline="") as f:
        summary_rows = list(csv.DictReader(f))
    assert len(summary_rows) == 1
    assert list(summary_rows[0].keys()) == COLLAPSE_SUMMARY_COLUMNS


def test_collapse_detection_uses_distance_and_consecutive_step_threshold(
    tmp_path: Path,
) -> None:
    cfg = _make_cfg(
        motor_torque_Nm=0.0, hook_enabled=True, n_flagella=2
    ).with_overrides(
        {
            "diagnostics": {
                "collapse": {
                    "enabled": True,
                    "write_every_step": True,
                    "collapse_distance_um": 0.15,
                    "collapse_consecutive_steps": 3,
                }
            }
        }
    )
    sim = Simulator(cfg)
    from sim_swim.sim.collapse_diagnostics import CollapseDiagnosticsRecorder

    rec = CollapseDiagnosticsRecorder(sim.model, cfg, tmp_path / "sim")
    pos_far = sim.model.positions_m.copy()
    rec.record(0, 0.0, _dummy_diag(pos_far))

    pos_close = pos_far.copy()
    f0 = sim.model.flagella_indices[0]
    f1 = sim.model.flagella_indices[1]
    pos_close[int(f1[0])] = pos_close[int(f0[0])] + np.array([1.0e-8, 0.0, 0.0])
    rec.record(1, 0.001, _dummy_diag(pos_close))
    rec.record(2, 0.002, _dummy_diag(pos_close))
    rec.record(3, 0.003, _dummy_diag(pos_close))
    _, summary_path = rec.write_files()

    with summary_path.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 1
    row = rows[0]
    assert row["collapse_detected"] in {"True", "true", "1"}
    assert int(row["first_collapse_step"]) == 1
    assert float(row["first_collapse_t_s"]) == pytest.approx(0.001)


def test_local_helix_is_not_used_when_template_projection_enabled(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = _make_cfg(
        motor_torque_Nm=0.0,
        hook_enabled=True,
        n_flagella=1,
    ).with_overrides(
        {
            "projection": {
                "enable_flagella_template_projection": True,
                "enable_flagella_chain_length_projection_when_template_off": True,
            },
            "local_helix": {"enabled": True},
        }
    )
    sim = Simulator(cfg)

    calls = {"count": 0}

    def _count_calls(positions_m: np.ndarray) -> np.ndarray:
        calls["count"] += 1
        return positions_m

    monkeypatch.setattr(sim.engine, "_project_local_helix_constraint", _count_calls)
    sim.engine._project_hook_and_flag_bonds(sim.model.positions_m.copy())
    assert calls["count"] == 0


def test_local_helix_is_used_when_template_projection_disabled(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    cfg = _make_cfg(
        motor_torque_Nm=0.0,
        hook_enabled=True,
        n_flagella=1,
    ).with_overrides(
        {
            "projection": {
                "enable_flagella_template_projection": False,
                "enable_flagella_chain_length_projection_when_template_off": True,
            },
            "local_helix": {"enabled": True},
        }
    )
    sim = Simulator(cfg)

    calls = {"count": 0}

    def _count_calls(positions_m: np.ndarray) -> np.ndarray:
        calls["count"] += 1
        return positions_m

    monkeypatch.setattr(sim.engine, "_project_local_helix_constraint", _count_calls)
    sim.engine._project_hook_and_flag_bonds(sim.model.positions_m.copy())
    assert calls["count"] > 0


def test_local_helix_targets_are_limited_to_n_local() -> None:
    cfg = _make_cfg(
        motor_torque_Nm=0.0,
        hook_enabled=True,
        n_flagella=1,
    ).with_overrides(
        {
            "local_helix": {
                "enabled": True,
                "n_local": 4,
            }
        }
    )
    sim = Simulator(cfg)

    assert len(sim.engine.local_helix_target_global_indices) == 1
    targets = sim.engine.local_helix_target_global_indices[0]
    assert int(targets.shape[0]) == 4
    expected = sim.model.flagella_indices[0].astype(int, copy=False)[:4]
    assert np.array_equal(targets, expected)
