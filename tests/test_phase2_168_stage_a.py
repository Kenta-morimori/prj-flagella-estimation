from __future__ import annotations

import csv
from pathlib import Path

import pytest
import yaml

from sim_swim.analysis.cli_profiles import args_from_profile, load_profile
from sim_swim.analysis.stage_a_2015_analysis import (
    THRESHOLD_POLICY,
    evaluate_motor_on,
    propose_thresholds,
)
from sim_swim.analysis.sweeps import stage_a_2015
from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig

ROOT = Path(__file__).parents[1]


def _load_config(name: str) -> dict:
    return yaml.safe_load((ROOT / "conf" / name).read_text(encoding="utf-8"))


@pytest.mark.parametrize(
    ("profile_name", "stage", "expected_steps"),
    [
        ("2015_stage_a_motor_off.yaml", "motor_off", 10_000),
        ("2015_stage_a_motor_on.yaml", "motor_on", 100_000),
    ],
)
def test_stage_a_profiles_fix_the_accepted_step_contract(
    profile_name: str,
    stage: str,
    expected_steps: int,
) -> None:
    profile = load_profile(ROOT / "conf/phase2_sweeps" / profile_name)
    args = stage_a_2015._parse_args(args_from_profile(profile))

    assert args.stage == stage
    assert args.profiles == ["project", "paper"]
    contract = stage_a_2015.STAGE_CONTRACT[stage]
    for config_name in ("sim_swim_2015.yaml", "sim_swim_2015_paper.yaml"):
        cfg = SimulationConfig.from_dict(_load_config(config_name)).with_overrides(
            {
                "time": {
                    "duration": {"value": contract["duration_tau"], "unit": "tau"},
                    "integration": {"dt_star": 1.0e-5},
                },
                "motor": {
                    "enabled": contract["motor_enabled"],
                    "enable_switching": False,
                },
            }
        )
        assert cfg.total_steps == expected_steps
        assert cfg.motor.enabled is contract["motor_enabled"]
        assert cfg.motor.reference_torque_Nm == pytest.approx(1.2e-18)


def test_simulator_can_sample_states_without_sampling_step_diagnostics(
    tmp_path: Path,
) -> None:
    cfg = SimulationConfig.from_dict(
        _load_config("sim_swim_2015_paper.yaml")
    ).with_overrides({"time": {"duration": {"value": 3.0e-5, "unit": "tau"}}})
    states = Simulator(cfg).run(
        cfg.time.duration_s,
        step_summary_dir=tmp_path,
        record_body_diagnostics=False,
        record_body_local_diagnostics=False,
        state_sample_interval_steps=2,
    )

    assert len(states) == 3
    assert states[0].t == pytest.approx(0.0)
    assert states[-1].t == pytest.approx(cfg.time.duration_s)
    with (tmp_path / "step_summary.csv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 3
    for field in (
        "flag_helix_radius_abs_err_over_b_max",
        "flag_helix_pitch_rel_err_max",
        "motor_force_balance_residual_ratio",
        "motor_torque_balance_residual_ratio",
        "net_force_residual_ratio",
        "net_torque_residual_ratio",
    ):
        assert field in rows[0]
        assert float(rows[0][field]) >= 0.0


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _write_analysis_fixture(root: Path, *, motor_on: bool) -> None:
    summary_rows = []
    for profile in ("project", "paper"):
        summary_rows.append(
            {
                "condition_id": profile,
                "profile": profile,
                "status": "completed",
                "completion_pass": True,
                "finite_pass_all": True,
                "axis_center_body_relative_net_abs_revolutions_max": 0.0001,
                "axis_center_body_relative_net_abs_revolutions_min": (
                    0.02 if motor_on else 0.0001
                ),
                "max_motor_force_balance_residual_ratio": 1.0e-12,
                "max_motor_torque_balance_residual_ratio": 1.0e-12,
                "body_shape_pass": True,
            }
        )
        _write_csv(
            root / profile / "step_summary.csv",
            [
                {
                    "flag_bond_rel_err_max": 0.001,
                    "hook_len_rel_err_max": 0.001,
                    "hook_angle_err_max_deg": 0.5,
                    "flag_bend_err_max_deg": 0.5,
                    "flag_torsion_err_max_deg": 0.5,
                    "flag_helix_radius_abs_err_over_b_max": 0.001,
                    "flag_helix_pitch_rel_err_max": 0.001,
                    "net_force_residual_ratio": 1.0e-12,
                    "net_torque_residual_ratio": 1.0e-12,
                    "motor_net_force_body_norm_N": 1.0e-20 if motor_on else 0.0,
                    "motor_net_force_flag_norm_N": 1.0e-20 if motor_on else 0.0,
                },
                {
                    "flag_bond_rel_err_max": 0.002,
                    "hook_len_rel_err_max": 0.002,
                    "hook_angle_err_max_deg": 1.0,
                    "flag_bend_err_max_deg": 1.0,
                    "flag_torsion_err_max_deg": 1.0,
                    "flag_helix_radius_abs_err_over_b_max": 0.002,
                    "flag_helix_pitch_rel_err_max": 0.002,
                    "net_force_residual_ratio": 2.0e-12,
                    "net_torque_residual_ratio": 2.0e-12,
                    "motor_net_force_body_norm_N": 1.0e-20 if motor_on else 0.0,
                    "motor_net_force_flag_norm_N": 1.0e-20 if motor_on else 0.0,
                },
            ],
        )
        _write_csv(
            root / profile / "body_constraint_diagnostics.csv",
            [
                {
                    "body_spring_max_stretch_ratio": 0.001,
                    "body_length_um": 2.0,
                    "body_width_mean_um": 1.0,
                    "body_width_min_um": 1.0,
                    "body_width_max_um": 1.0,
                    "body_cross_section_area_min_um2": 1.0,
                    "body_cross_section_area_max_um2": 1.0,
                },
                {
                    "body_spring_max_stretch_ratio": 0.002,
                    "body_length_um": 2.002,
                    "body_width_mean_um": 1.001,
                    "body_width_min_um": 0.999,
                    "body_width_max_um": 1.001,
                    "body_cross_section_area_min_um2": 0.999,
                    "body_cross_section_area_max_um2": 1.001,
                },
            ],
        )
    _write_csv(root / "summary.csv", summary_rows)


def test_motor_off_analysis_proposes_shared_thresholds(tmp_path: Path) -> None:
    motor_off = tmp_path / "motor_off"
    _write_analysis_fixture(motor_off, motor_on=False)

    proposal = propose_thresholds(motor_off)

    assert proposal["status"] == "proposed"
    assert proposal["failures"] == []
    assert proposal["thresholds"]["max_flag_bond_rel_err"] == pytest.approx(0.01)
    assert proposal["thresholds"][
        "min_axis_center_body_relative_revolutions"
    ] == pytest.approx(0.01)


def test_motor_on_analysis_requires_rotation_and_visual_review(tmp_path: Path) -> None:
    motor_on = tmp_path / "motor_on"
    _write_analysis_fixture(motor_on, motor_on=True)
    thresholds = {metric: policy["cap"] for metric, policy in THRESHOLD_POLICY.items()}
    thresholds.update(
        {
            "min_axis_center_body_relative_revolutions": 0.01,
            "max_motor_force_balance_residual_ratio": 1.0e-8,
            "max_motor_torque_balance_residual_ratio": 1.0e-8,
        }
    )

    decision = evaluate_motor_on(motor_on, {"thresholds": thresholds})

    assert decision["next_action"] == "perform_visual_review"
    assert all(item["automated_pass"] for item in decision["profiles"].values())
    assert all(item["visual_review_required"] for item in decision["profiles"].values())
