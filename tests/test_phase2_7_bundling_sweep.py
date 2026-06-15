from __future__ import annotations

import importlib.util
from pathlib import Path


SCRIPT_PATH = (
    Path(__file__).resolve().parents[1]
    / "scripts"
    / "01_simulate_swimming"
    / "run_phase2_7_bundling_sweep.py"
)
SPEC = importlib.util.spec_from_file_location("phase27_bundling_sweep", SCRIPT_PATH)
assert SPEC is not None
phase27 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(phase27)


def test_phase27_classifies_shape_fail_as_collapse() -> None:
    row = {"shape_pass_nonbody": "False"}

    assert phase27.classify_phase27_condition(row, n_flagella=3) == "collapse"


def test_phase27_classifies_posterior_bundle() -> None:
    row = {
        "shape_pass_nonbody": "True",
        "bundle_axis_vs_rear_angle_deg": "12.0",
        "bundle_rearward_projection": "0.9",
        "bundle_participation_ratio": "1.0",
        "bundle_independent_flagella_count": "0",
    }

    assert phase27.classify_phase27_condition(row, n_flagella=3) == "posterior_bundle"


def test_phase27_classifies_partial_bundle() -> None:
    row = {
        "shape_pass_nonbody": "True",
        "bundle_axis_vs_rear_angle_deg": "20.0",
        "bundle_rearward_projection": "0.8",
        "bundle_participation_ratio": "0.67",
        "bundle_independent_flagella_count": "1",
    }

    assert phase27.classify_phase27_condition(row, n_flagella=3) == "partial_bundle"


def test_phase27_classifies_nonposterior_condition_as_no_bundle() -> None:
    row = {
        "shape_pass_nonbody": "True",
        "bundle_axis_vs_rear_angle_deg": "90.0",
        "bundle_rearward_projection": "0.0",
        "bundle_participation_ratio": "1.0",
        "bundle_independent_flagella_count": "0",
    }

    assert phase27.classify_phase27_condition(row, n_flagella=3) == "no_bundle"


def test_phase27_net_abs_helix_spin_revolutions_uses_phase_delta() -> None:
    rows = [
        {"flag_helix_spin_phase_deg": "10.0"},
        {"flag_helix_spin_phase_deg": "370.0"},
    ]

    assert phase27._net_abs_helix_spin_revolutions(rows) == 1.0


def test_phase27_build_config_sets_initial_tangent_angle() -> None:
    raw_cfg = phase27._load_yaml(
        Path(__file__).resolve().parents[1] / "conf/sim_swim.yaml"
    )

    cfg = phase27._build_config(
        raw_cfg,
        orientation_mode="side_attach",
        initial_flagellum_axis_from_rear_deg=10.0,
        n_flagella=3,
        torque_Nm=0.5e-20,
        duration_s=0.02,
        dt_star=1.0e-4,
        overrides=[],
    )

    assert cfg.flagella.initial_orientation_mode == "side_attach"
    assert cfg.flagella.initial_flagellum_axis_from_rear_deg == 10.0
    assert cfg.flagella.n_flagella == 3
    assert cfg.motor.torque_Nm == 0.5e-20
