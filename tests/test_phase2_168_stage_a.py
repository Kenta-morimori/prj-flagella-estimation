from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
import pytest
import yaml

from sim_swim.analysis.flagella_count_behavior import save_state_archive
from sim_swim.analysis.cli_profiles import args_from_profile, load_profile
from sim_swim.analysis import stage_a_2015_analysis
from sim_swim.analysis.stage_a_2015_analysis import (
    THRESHOLD_POLICY,
    evaluate_motor_on,
    propose_thresholds,
)
from sim_swim.analysis.sweeps import stage_a_2015
from sim_swim.sim.core import Simulator
from sim_swim.sim.core import SimulationState
from sim_swim.sim.params import SimulationConfig

ROOT = Path(__file__).parents[1]


def _load_config(name: str) -> dict:
    return yaml.safe_load((ROOT / "conf" / name).read_text(encoding="utf-8"))


@pytest.mark.parametrize(
    ("profile_name", "stage", "expected_steps", "dt_star", "duration_tau"),
    [
        ("2015_stage_a_motor_off.yaml", "motor_off", 10_000, 1.0e-5, 0.1),
        ("2015_stage_a_motor_on.yaml", "motor_on", 100_000, 1.0e-5, 1.0),
        (
            "2015_stage_a_dt1e4_motor_off_reference.yaml",
            "motor_off",
            1_000,
            1.0e-4,
            0.1,
        ),
        (
            "2015_stage_a_dt1e4_motor_on_reference.yaml",
            "motor_on",
            10_000,
            1.0e-4,
            1.0,
        ),
        (
            "2015_stage_a_dt1e5_motor_on_short_reference.yaml",
            "motor_on",
            10_000,
            1.0e-5,
            0.1,
        ),
    ],
)
def test_stage_a_profiles_fix_the_accepted_step_contract(
    profile_name: str,
    stage: str,
    expected_steps: int,
    dt_star: float,
    duration_tau: float,
) -> None:
    profile = load_profile(ROOT / "conf/phase2_sweeps" / profile_name)
    args = stage_a_2015._parse_args(args_from_profile(profile))

    assert args.stage == stage
    assert args.profiles == ["project", "paper"]
    contract = stage_a_2015.STAGE_CONTRACT[stage]
    effective_duration = (
        contract["duration_tau"] if args.duration_tau is None else args.duration_tau
    )
    assert args.dt_star == pytest.approx(dt_star)
    assert effective_duration == pytest.approx(duration_tau)
    for config_name in ("sim_swim_2015.yaml", "sim_swim_2015_paper.yaml"):
        cfg = SimulationConfig.from_dict(_load_config(config_name)).with_overrides(
            {
                "time": {
                    "duration": {"value": effective_duration, "unit": "tau"},
                    "integration": {"dt_star": args.dt_star},
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


def test_stage_a_smoke_duration_uses_selected_dt_star(
    capsys: pytest.CaptureFixture[str],
) -> None:
    stage_a_2015.run_stage_a(
        [
            "--stage",
            "motor_on",
            "--dt-star",
            "1e-4",
            "--smoke-steps",
            "3",
            "--dry-run",
        ]
    )

    assert "duration_tau=0.00030000000000000003" in capsys.readouterr().out


def test_torque_screen_profile_covers_two_orders_at_2000_steps(
    capsys: pytest.CaptureFixture[str],
) -> None:
    profile = load_profile(ROOT / "conf/phase2_sweeps/2015_stage_a_torque_screen.yaml")
    args = stage_a_2015._parse_args(args_from_profile(profile))

    assert args.duration_tau == pytest.approx(0.02)
    assert args.dt_star == pytest.approx(1.0e-5)
    assert args.motor_torque_scales == [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0]
    assert int(round(args.duration_tau / args.dt_star)) == 2_000

    stage_a_2015.run_stage_a(args_from_profile(profile) + ["--dry-run"])
    output = capsys.readouterr().out
    assert output.count("torque_screen_dt1e5_short") == 14
    assert "project_torque_x0p1" in output
    assert "paper_torque_x10" in output


def test_torque_screen_condition_ids_are_unique() -> None:
    ids = {
        stage_a_2015._condition_id(profile, scale, 7)
        for profile in ("project", "paper")
        for scale in (0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0)
    }

    assert len(ids) == 14


def test_project_torque_long_profile_has_three_1tau_conditions(
    capsys: pytest.CaptureFixture[str],
) -> None:
    profile = load_profile(
        ROOT / "conf/phase2_sweeps/2015_stage_a_project_torque_long.yaml"
    )
    args = stage_a_2015._parse_args(args_from_profile(profile))

    assert args.profiles == ["project"]
    assert args.duration_tau == pytest.approx(1.0)
    assert args.dt_star == pytest.approx(1.0e-5)
    assert args.motor_torque_scales == [0.5, 1.0, 2.0]

    stage_a_2015.run_stage_a(args_from_profile(profile) + ["--dry-run"])
    output = capsys.readouterr().out
    assert output.count("project_torque_long_dt1e5") == 3
    assert "project_torque_x0p5" in output
    assert "project_torque_x2" in output


def test_project_torque_dt1e4_grid_profile_covers_only_missing_cells(
    capsys: pytest.CaptureFixture[str],
) -> None:
    profile = load_profile(
        ROOT / "conf/phase2_sweeps/2015_stage_a_project_torque_dt1e4_missing.yaml"
    )
    args = stage_a_2015._parse_args(args_from_profile(profile))

    assert args.profiles == ["project"]
    assert args.duration_tau == pytest.approx(1.0)
    assert args.dt_star == pytest.approx(1.0e-4)
    assert args.motor_torque_scales == [0.5, 2.0]
    project_cfg = SimulationConfig.from_dict(_load_config("sim_swim_2015.yaml"))
    assert project_cfg.flagella.initial_helix_axis_from_rear_deg == pytest.approx(0.0)

    stage_a_2015.run_stage_a(args_from_profile(profile) + ["--dry-run"])
    output = capsys.readouterr().out
    assert output.count("project_torque_dt1e4_grid_missing") == 2
    assert "project_torque_x0p5" in output
    assert "project_torque_x2" in output
    assert "project_torque_x1" not in output


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
    assert "min_axis_center_body_relative_revolutions" not in proposal["thresholds"]


def test_motor_on_analysis_requires_rotation_and_visual_review(tmp_path: Path) -> None:
    motor_on = tmp_path / "motor_on"
    _write_analysis_fixture(motor_on, motor_on=True)
    thresholds = {metric: policy["cap"] for metric, policy in THRESHOLD_POLICY.items()}
    thresholds.update(
        {
            "max_motor_force_balance_residual_ratio": 1.0e-8,
            "max_motor_torque_balance_residual_ratio": 1.0e-8,
        }
    )

    decision = evaluate_motor_on(motor_on, {"thresholds": thresholds})

    assert decision["next_action"] == "perform_visual_review"
    assert all(item["automated_pass"] for item in decision["profiles"].values())
    assert all(item["visual_review_required"] for item in decision["profiles"].values())


def _state(t: float, *, shift: float = 0.0) -> SimulationState:
    positions = np.asarray(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]], dtype=float
    )
    positions[:, 1] += shift
    return SimulationState(
        t=t,
        position_um=(0.0, 0.0, 0.0),
        quaternion=(0.0, 0.0, 0.0, 1.0),
        velocity_um_s=(0.0, 0.0, 0.0),
        omega_rad_s=(0.0, 0.0, 0.0),
        bead_positions_um=positions,
        flag_states=(0,),
        reverse_flagella=(),
    )


def _write_axis_rows(path: Path, final_phase: float) -> None:
    _write_csv(
        path,
        [
            {
                "t_s": 0.1 * fraction,
                "flag_id": 0,
                "axis_center_body_relative_phase_deg": final_phase * fraction,
            }
            for fraction in (0.0, 1.0 / 3.0, 2.0 / 3.0, 1.0)
        ],
    )


def test_paired_dt_comparison_passes_and_detects_rotation_reversal(
    tmp_path: Path,
) -> None:
    coarse = tmp_path / "coarse"
    fine = tmp_path / "fine"
    for root in (coarse, fine):
        save_state_archive(
            root / "project/state_archive.npz", [_state(0.0), _state(0.1)]
        )
    _write_axis_rows(
        coarse / "project/flag_helix_axis_diagnostics.csv", final_phase=360.0
    )
    _write_axis_rows(
        fine / "project/flag_helix_axis_diagnostics.csv", final_phase=350.0
    )
    manifest = {
        "conditions": [
            {
                "condition_id": "project",
                "comparison_scales": {"b_um": 1.0, "body_beads": 2},
            }
        ]
    }

    passed = stage_a_2015_analysis._paired_profile_comparison(
        profile="project",
        coarse_run=coarse,
        fine_run=fine,
        coarse_manifest=manifest,
        fine_manifest=manifest,
    )
    assert passed["paired_pass"] is True

    _write_axis_rows(
        coarse / "project/flag_helix_axis_diagnostics.csv", final_phase=-360.0
    )
    failed = stage_a_2015_analysis._paired_profile_comparison(
        profile="project",
        coarse_run=coarse,
        fine_run=fine,
        coarse_manifest=manifest,
        fine_manifest=manifest,
    )
    assert failed["paired_pass"] is False
    assert failed["flag_rotation"][0]["same_direction"] is False
