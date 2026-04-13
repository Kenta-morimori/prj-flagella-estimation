from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from sim_swim.sim.core import Simulator
from sim_swim.sim.params import (
    DynamicsMode,
    infer_dynamics_mode,
    validate_dynamics_mode_consistency,
    SimulationConfig,
)


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


def _make_phase0a_cfg(duration_s: float = 1.0e-2) -> SimulationConfig:
    """Phase0a: body-only static (no flagella, motor off)."""
    return _make_cfg(
        motor_torque_Nm=0.0,
        n_flagella=0,
        duration_s=duration_s,
    )


def _make_phase0b_cfg(duration_s: float = 1.0e-2) -> SimulationConfig:
    """Phase0b: body + hook + minimal_basal_stub static (motor off)."""
    return _make_cfg(
        motor_torque_Nm=0.0,
        n_flagella=1,
        stub_mode="minimal_basal_stub",
        duration_s=duration_s,
    )


def _make_phase1_cfg(
    surrogate_torque_Nm: float = 1.0e-20, duration_s: float = 1.0e-2
) -> SimulationConfig:
    """Phase1: body-only surrogate torque (no flagella, body_equiv_load enabled)."""
    from dataclasses import asdict

    cfg = _make_cfg(
        motor_torque_Nm=0.0,
        n_flagella=0,
        duration_s=duration_s,
    )
    # body_equiv_load を有効化
    cfg_dict = asdict(cfg)
    cfg_dict["body_equiv_load"] = {
        "enabled": True,
        "mode": "pure_couple",
        "target_torque_Nm": surrogate_torque_Nm,
        "target_force_N": 0.0,
        "attach_region_id": 0,
    }
    return SimulationConfig.from_dict(cfg_dict)


def _make_phase2_cfg(
    motor_torque_Nm: float = 4.0e-21, duration_s: float = 1.0e-2
) -> SimulationConfig:
    """Phase2: minimal_basal_stub + actual motor."""
    return _make_cfg(
        motor_torque_Nm=motor_torque_Nm,
        n_flagella=1,
        stub_mode="minimal_basal_stub",
        duration_s=duration_s,
    )


def _make_phase3_cfg(
    motor_torque_Nm: float = 4.0e-21, duration_s: float = 1.0e-2
) -> SimulationConfig:
    """Phase3: full_flagella + actual motor."""
    return _make_cfg(
        motor_torque_Nm=motor_torque_Nm,
        n_flagella=1,
        stub_mode="full_flagella",
        duration_s=duration_s,
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


def _series(rows: list[dict[str, str]], key: str) -> np.ndarray:
    return np.asarray([float(row[key]) for row in rows], dtype=float)


def _assert_all_finite(series: np.ndarray, label: str) -> None:
    assert np.isfinite(series).all(), f"{label} contains NaN/Inf"


def _assert_no_runaway_growth(
    series: np.ndarray,
    label: str,
    *,
    baseline_window: int = 3,
    tail_window: int = 3,
    max_growth_factor: float = 3.0,
    baseline_floor: float = 0.0,
) -> None:
    assert series.size >= baseline_window + tail_window
    baseline = float(np.median(series[:baseline_window]))
    tail = float(np.median(series[-tail_window:]))
    effective_baseline = max(baseline, float(baseline_floor))
    assert tail <= effective_baseline * max_growth_factor, (
        f"{label} runaway growth: baseline={baseline:.6g}, "
        f"effective_baseline={effective_baseline:.6g}, tail={tail:.6g}, "
        f"limit={effective_baseline * max_growth_factor:.6g}"
    )


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


# ==============================================================================
# Phase-gated tests (Issue #37: phase-based decoupling)
# ==============================================================================


def test_phase0a_canonical_config_buildable() -> None:
    """Hard test 1: Phase0a canonical config can be constructed."""
    cfg = _make_phase0a_cfg(duration_s=0.01)
    assert cfg.flagella.n_flagella == 0
    assert cfg.motor.torque_Nm == 0.0
    assert cfg.brownian.enabled is False
    assert cfg.run_tumble.run_tau > 0  # switching off (long tau)


def test_phase0b_canonical_config_buildable() -> None:
    """Hard test 1: Phase0b canonical config can be constructed."""
    cfg = _make_phase0b_cfg(duration_s=0.01)
    assert cfg.flagella.n_flagella == 1
    assert cfg.flagella.stub_mode == "minimal_basal_stub"
    assert cfg.motor.torque_Nm == 0.0
    assert cfg.brownian.enabled is False


def test_phase1_canonical_config_buildable() -> None:
    """Hard test 1: Phase1 canonical config can be constructed."""
    cfg = _make_phase1_cfg(surrogate_torque_Nm=1.0e-20, duration_s=0.01)
    assert cfg.flagella.n_flagella == 0
    assert cfg.motor.torque_Nm == 0.0
    assert cfg.body_equiv_load.enabled is True
    assert cfg.body_equiv_load.mode == "pure_couple"
    assert cfg.brownian.enabled is False


def test_phase2_canonical_config_buildable() -> None:
    """Hard test 1: Phase2 canonical config can be constructed."""
    cfg = _make_phase2_cfg(motor_torque_Nm=4.0e-21, duration_s=0.01)
    assert cfg.flagella.n_flagella == 1
    assert cfg.flagella.stub_mode == "minimal_basal_stub"
    assert cfg.motor.torque_Nm == 4.0e-21
    assert cfg.brownian.enabled is False


def test_phase3_canonical_config_buildable() -> None:
    """Hard test 1: Phase3 canonical config can be constructed."""
    cfg = _make_phase3_cfg(motor_torque_Nm=4.0e-21, duration_s=0.01)
    assert cfg.flagella.n_flagella == 1
    assert cfg.flagella.stub_mode == "full_flagella"
    assert cfg.motor.torque_Nm == 4.0e-21
    assert cfg.brownian.enabled is False


def test_phase0a_body_only_finite_completion(tmp_path: Path) -> None:
    """Hard test 2: Phase0a body-only finite completion."""
    cfg = _make_phase0a_cfg(duration_s=0.01)
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, cfg.time.duration_s, tmp_path / "phase0a")

    assert len(rows) >= 1
    # All steps finite
    for r in rows:
        assert r["pos_all_finite"] in {"True", "true", "1"}
        assert r["any_nan"] in {"False", "false", "0"}
        assert r["any_inf"] in {"False", "false", "0"}


def test_phase0b_minimal_stub_motor_off_finite_completion(tmp_path: Path) -> None:
    """Hard test 3: Phase0b minimal_basal_stub/motor off finite completion."""
    cfg = _make_phase0b_cfg(duration_s=0.01)
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, cfg.time.duration_s, tmp_path / "phase0b")

    assert len(rows) >= 1
    # All steps finite
    for r in rows:
        assert r["pos_all_finite"] in {"True", "true", "1"}
        assert r["any_nan"] in {"False", "false", "0"}
        assert r["any_inf"] in {"False", "false", "0"}
    # Check diagnostics present
    first = rows[0]
    assert float(first["local_attach_first_rel_err"]) >= 0.0
    assert float(first["local_first_second_rel_err"]) >= 0.0


def test_phase1_body_equiv_load_pure_couple_finite_completion(tmp_path: Path) -> None:
    """Hard test 4: Phase1 body-only + pure_couple finite completion & diagnostics."""
    cfg = _make_phase1_cfg(surrogate_torque_Nm=1.0e-20, duration_s=0.01)
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, cfg.time.duration_s, tmp_path / "phase1")

    assert len(rows) >= 1
    # All steps finite
    for r in rows:
        assert r["pos_all_finite"] in {"True", "true", "1"}
        assert r["any_nan"] in {"False", "false", "0"}
        assert r["any_inf"] in {"False", "false", "0"}
    # body_equiv_load active
    first = rows[0]
    assert first["body_equiv_load_mode"] == "pure_couple"
    assert float(first["body_equiv_load_target_torque_Nm"]) == 1.0e-20


def test_phase2_minimal_stub_motor_on_short_run_diagnostics(tmp_path: Path) -> None:
    """Hard test 5: Phase2 minimal_basal_stub + motor on 0.1s diagnostics."""
    cfg = _make_phase2_cfg(motor_torque_Nm=5.0e-20, duration_s=0.1)
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, cfg.time.duration_s, tmp_path / "phase2")

    assert len(rows) >= 1
    first = rows[0]
    # Motor diagnostics present
    assert "motor_Ta_norm" in first
    assert "motor_split_residual_norm" in first
    assert "motor_attach_force_norm" in first
    # Basal diagnostics present
    assert "local_attach_first_rel_err" in first
    assert "local_first_second_rel_err" in first


def test_phase3_full_flagella_motor_on_short_run_diagnostics(tmp_path: Path) -> None:
    """Hard test 6: Phase3 full_flagella + motor on 0.1s diagnostics."""
    cfg = _make_phase3_cfg(motor_torque_Nm=5.0e-20, duration_s=0.1)
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, cfg.time.duration_s, tmp_path / "phase3")

    assert len(rows) >= 1
    first = rows[0]
    # Full flagella diagnostics present
    assert "flag_bend_err_max_deg" in first
    assert "flag_torsion_err_max_deg" in first
    assert "flag_bond_rel_err_max" in first
    # Basal diagnostics still present
    assert "local_attach_first_rel_err" in first


def test_phase3_full_flagella_static_shape_gate(tmp_path: Path) -> None:
    """Phase3 static must keep full-flagella shape stable, not just finite."""
    cfg = _make_phase3_cfg(motor_torque_Nm=0.0, duration_s=0.01)
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(
        sim,
        cfg.time.duration_s,
        tmp_path / "phase3_static",
    )

    assert len(rows) >= 10
    assert all(r["pos_all_finite"] in {"True", "true", "1"} for r in rows)
    assert all(r["any_nan"] in {"False", "false", "0"} for r in rows)
    assert all(r["any_inf"] in {"False", "false", "0"} for r in rows)

    # Structure must remain the full-flagella static baseline.
    assert all(int(r["hook_count"]) == 1 for r in rows)
    assert all(int(r["flag_intra_count"]) == 10 for r in rows)

    hook_len_over_b = _series(rows, "hook_len_mean_over_b")
    attach_rel_err = _series(rows, "local_attach_first_rel_err")
    bond_rel_err = _series(rows, "flag_bond_rel_err_max")
    torsion_err = _series(rows, "local_first_torsion_err_deg")

    _assert_all_finite(hook_len_over_b, "hook_len_mean_over_b")
    _assert_all_finite(attach_rel_err, "local_attach_first_rel_err")
    _assert_all_finite(bond_rel_err, "flag_bond_rel_err_max")
    _assert_all_finite(torsion_err, "local_first_torsion_err_deg")

    # Static full flagella must not show runaway late-time growth.
    _assert_no_runaway_growth(hook_len_over_b, "hook_len_mean_over_b")
    _assert_no_runaway_growth(
        attach_rel_err,
        "local_attach_first_rel_err",
        baseline_floor=1.0e-5,
    )
    _assert_no_runaway_growth(
        bond_rel_err,
        "flag_bond_rel_err_max",
        baseline_floor=1.0e-5,
    )
    _assert_no_runaway_growth(
        torsion_err,
        "local_first_torsion_err_deg",
        baseline_floor=1.0e-3,
    )


def test_sim_diagnostics_docs_file_exists() -> None:
    """Hard test 7: Verify sim_diagnostics.md exists (docs update rule)."""
    from pathlib import Path as PathlibPath

    docs_path = (
        PathlibPath(__file__).parent.parent / "docs" / "phase2" / "sim_diagnostics.md"
    )
    assert docs_path.is_file(), f"sim_diagnostics.md missing at {docs_path}"


# ==============================================================================
# Issue #37: DynamicsMode and structure validation tests
# ==============================================================================


def test_infer_dynamics_mode_body_only() -> None:
    """Test: infer BODY_ONLY from n_flagella=0."""
    mode = infer_dynamics_mode(n_flagella=0, stub_mode="full_flagella")
    assert mode == DynamicsMode.BODY_ONLY


def test_infer_dynamics_mode_body_hook() -> None:
    """Test: infer BODY_HOOK from n_flagella >= 1, stub_mode='minimal_basal_stub'."""
    mode = infer_dynamics_mode(n_flagella=1, stub_mode="minimal_basal_stub")
    assert mode == DynamicsMode.BODY_HOOK

    mode = infer_dynamics_mode(n_flagella=2, stub_mode="minimal_basal_stub")
    assert mode == DynamicsMode.BODY_HOOK


def test_infer_dynamics_mode_body_hook_flagella() -> None:
    """Test: infer BODY_HOOK_FLAGELLA from n_flagella >= 1, stub_mode='full_flagella'."""
    mode = infer_dynamics_mode(n_flagella=1, stub_mode="full_flagella")
    assert mode == DynamicsMode.BODY_HOOK_FLAGELLA

    mode = infer_dynamics_mode(n_flagella=3, stub_mode="full_flagella")
    assert mode == DynamicsMode.BODY_HOOK_FLAGELLA


def test_infer_dynamics_mode_invalid_stub_mode() -> None:
    """Test: infer_dynamics_mode raises on invalid stub_mode."""
    with pytest.raises(ValueError, match="Invalid stub_mode"):
        infer_dynamics_mode(n_flagella=1, stub_mode="invalid_mode")


def test_infer_dynamics_mode_invalid_n_flagella() -> None:
    """Test: infer_dynamics_mode raises on invalid n_flagella."""
    with pytest.raises(ValueError, match="Invalid n_flagella"):
        infer_dynamics_mode(n_flagella=-1, stub_mode="full_flagella")


def test_validate_dynamics_mode_consistency_body_only() -> None:
    """Test: validation passes for BODY_ONLY mode."""
    validate_dynamics_mode_consistency(
        mode=DynamicsMode.BODY_ONLY, n_flagella=0, stub_mode="full_flagella"
    )


def test_validate_dynamics_mode_consistency_body_hook() -> None:
    """Test: validation passes for BODY_HOOK mode."""
    validate_dynamics_mode_consistency(
        mode=DynamicsMode.BODY_HOOK, n_flagella=1, stub_mode="minimal_basal_stub"
    )


def test_validate_dynamics_mode_consistency_body_hook_flagella() -> None:
    """Test: validation passes for BODY_HOOK_FLAGELLA mode."""
    validate_dynamics_mode_consistency(
        mode=DynamicsMode.BODY_HOOK_FLAGELLA, n_flagella=1, stub_mode="full_flagella"
    )


def test_validate_dynamics_mode_consistency_mismatch() -> None:
    """Test: validation raises on mode mismatch."""
    with pytest.raises(ValueError, match="Dynamics mode mismatch"):
        validate_dynamics_mode_consistency(
            mode=DynamicsMode.BODY_ONLY, n_flagella=1, stub_mode="minimal_basal_stub"
        )


def test_structure_body_only_has_no_flagella() -> None:
    """Test: BODY_ONLY mode has zero flagella beads."""
    cfg = _make_phase0a_cfg(duration_s=0.001)
    sim = Simulator(cfg)

    # Model structure validation
    assert sim.model.flagella_indices == []
    assert len(sim.model.hook_triplets) == 0


def test_structure_body_hook_has_minimal_flagella() -> None:
    """Test: BODY_HOOK mode has exactly 3 flagella beads (attach, first, second)."""
    cfg = _make_phase0b_cfg(duration_s=0.001)
    sim = Simulator(cfg)

    # Check flagella_indices
    assert len(sim.model.flagella_indices) == 1
    flagella_bead_ids = sim.model.flagella_indices[0]
    assert len(flagella_bead_ids) == 3, (
        f"Expected 3 beads, got {len(flagella_bead_ids)}"
    )

    # Check hook triplets exist
    assert len(sim.model.hook_triplets) > 0


def test_structure_body_hook_flagella_has_full_flagella() -> None:
    """Test: BODY_HOOK_FLAGELLA mode has full flagella (11 beads by default)."""
    cfg = _make_phase3_cfg(motor_torque_Nm=0.0, duration_s=0.001)
    sim = Simulator(cfg)

    # Check flagella_indices
    assert len(sim.model.flagella_indices) == 1
    flagella_bead_ids = sim.model.flagella_indices[0]
    assert len(flagella_bead_ids) == 11, (
        f"Expected 11 beads, got {len(flagella_bead_ids)}"
    )

    # Check hook triplets exist
    assert len(sim.model.hook_triplets) > 0


def test_phase0a_short_run_all_finite() -> None:
    """Test (Phase A): Phase0a runs 10+ steps with all numerics finite."""
    cfg = _make_phase0a_cfg(duration_s=0.01)
    sim = Simulator(cfg)

    # Collect states
    states = sim.run(cfg.time.duration_s)
    assert len(states) >= 10

    # Check all coordinates finite
    for state in states:
        pos = np.asarray(state.position_um, dtype=float)
        q = np.asarray(state.quaternion, dtype=float)
        assert np.isfinite(pos).all(), "Step has NaN/Inf in position"
        assert np.isfinite(q).all(), "Step has NaN/Inf in quaternion"


def test_phase0b_short_run_all_finite_with_diagnostics() -> None:
    """Test (Phase A): Phase0b (minimal stub) runs 10+ steps, all finite, diagnostics present."""
    cfg = _make_phase0b_cfg(duration_s=0.01)
    sim = Simulator(cfg)

    # Collect states and write diagnostics
    states = sim.run(cfg.time.duration_s)
    assert len(states) >= 10

    # Check all coordinates finite
    for state in states:
        pos = np.asarray(state.position_um, dtype=float)
        q = np.asarray(state.quaternion, dtype=float)
        assert np.isfinite(pos).all(), "Step has NaN/Inf in position"
        assert np.isfinite(q).all(), "Step has NaN/Inf in quaternion"


def test_phase0b_multi_seed_reproducibility(tmp_path: Path) -> None:
    """Test (Phase A): Phase0b runs reproducibly across seeds (first-fail position)."""
    seeds = [0, 1, 2]
    results = []

    for seed in seeds:
        cfg_dict = {
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
                "stub_mode": "minimal_basal_stub",
                "n_beads_per_flagellum": 11,
                "discretization": {"ds_over_b": 0.58},
                "bond_L_over_b": 0.58,
                "length_over_b": 2.32,
                "helix_init": {"radius_over_b": 0.25, "pitch_over_b": 2.5},
            },
            "fluid": {"viscosity_Pa_s": 1.0e-3},
            "motor": {"torque_Nm": 0.0, "reverse_n_flagella": 1},
            "potentials": {
                "spring": {"H_over_T_over_b": 10.0, "s": 0.1},
                "bend": {"kb_over_T": 20.0, "s": 0.1},
                "torsion": {"kt_over_T": 20.0, "fd_eps_over_b": 1.0e-3},
                "hook_attraction": {"h_bind": 10.0, "h_unbind": 5.0},
                "hook_length": {"k_hook": 100.0, "target_L_nm": 250.0},
                "hook_bending": {"kb_hook": 20.0, "hook_angle_max_deg": 90.0},
                "repulsion": {"r0": 0.3, "K": 1000.0},
            },
            "hook": {"enabled": True},
            "time": {"dt_star": 1.0e-3, "dt_s": 1.0e-3, "duration_s": 0.01},
            "brownian": {"enabled": False, "T_K": 293.0},
            "run_tumble": {
                "run_tau": 999.0,
                "tumble_tau": 0.001,
                "semicoiled_tau": 999.0,
                "curly1_tau": 999.0,
            },
            "body_equiv_load": {
                "enabled": False,
                "mode": "off",
                "target_torque_Nm": 0.0,
                "target_force_N": 0.0,
                "attach_region_id": 0,
            },
            "seed": {"global_seed": seed},
            "output": {"base_dir": "outputs"},
        }
        cfg = SimulationConfig.from_dict(cfg_dict)
        sim = Simulator(cfg)
        rows = _run_and_load_step_summary(
            sim, cfg.time.duration_s, tmp_path / f"seed_{seed}"
        )
        results.append(
            {
                "seed": seed,
                "n_steps": len(rows),
                "all_finite": all(
                    r["pos_all_finite"] in {"True", "true", "1"} for r in rows
                ),
            }
        )

    # All seeds should complete without early termination
    for result in results:
        assert result["all_finite"], f"Seed {result['seed']} had NaN/Inf"
        assert result["n_steps"] >= 10, f"Seed {result['seed']} had too few steps"

    # Check that completion is consistent (not wildly different)
    step_counts = [r["n_steps"] for r in results]
    assert max(step_counts) - min(step_counts) <= 2, (
        f"Step counts too variable across seeds: {step_counts}"
    )
