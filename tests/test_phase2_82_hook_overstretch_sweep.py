from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

import pytest

from sim_swim.analysis.sweeps import hook_overstretch as script


def test_phase2_82_body_first_grid_conditions() -> None:
    args = SimpleNamespace(
        mode="body-first-grid",
        attach_first_spring_scales=[1.0, 2.0],
        body_axis_angle_scales=[1.0, 3.0],
        fixed_first_second_spring_scale=1.0,
    )

    conditions = script.build_conditions(args)

    assert [condition.condition_id for condition in conditions] == [
        "af1_axis1_fs1",
        "af1_axis3_fs1",
        "af2_axis1_fs1",
        "af2_axis3_fs1",
    ]
    assert all(condition.mode == "body-first-grid" for condition in conditions)
    assert conditions[-1].scales == {
        "local_attach_first_spring_scale": 2.0,
        "local_attach_first_body_axis_angle_scale": 3.0,
        "local_first_second_spring_scale": 1.0,
    }


def test_phase2_82_first_second_grid_conditions() -> None:
    args = SimpleNamespace(
        mode="first-second-grid",
        fixed_attach_first_spring_scale=2.0,
        fixed_body_axis_angle_scale=1.5,
        first_second_spring_scales=[1.0, 2.0, 3.0],
    )

    conditions = script.build_conditions(args)

    assert [condition.condition_id for condition in conditions] == [
        "af2_axis1p5_fs1",
        "af2_axis1p5_fs2",
        "af2_axis1p5_fs3",
    ]
    assert all(condition.mode == "first-second-grid" for condition in conditions)
    assert conditions[1].scales == {
        "local_attach_first_spring_scale": 2.0,
        "local_attach_first_body_axis_angle_scale": 1.5,
        "local_first_second_spring_scale": 2.0,
    }


def test_phase2_82_first_second_grid_can_fix_attach_frame_scales() -> None:
    args = SimpleNamespace(
        mode="first-second-grid",
        fixed_attach_first_spring_scale=1.0,
        fixed_body_axis_angle_scale=1.0,
        fixed_attach_frame_position_scale=3.0,
        fixed_attach_frame_tangent_scale=1.5,
        first_second_spring_scales=[1.0, 1.5],
    )

    conditions = script.build_conditions(args)

    assert [condition.condition_id for condition in conditions] == [
        "af1_axis1_fs1_fp3_ft1p5",
        "af1_axis1_fs1p5_fp3_ft1p5",
    ]
    assert conditions[-1].scales == {
        "local_attach_first_spring_scale": 1.0,
        "local_attach_first_body_axis_angle_scale": 1.0,
        "local_first_second_spring_scale": 1.5,
        "local_attach_frame_position_scale": 3.0,
        "local_attach_frame_tangent_scale": 1.5,
    }


def test_phase2_82_attach_frame_grid_conditions() -> None:
    args = SimpleNamespace(
        mode="attach-frame-grid",
        fixed_attach_first_spring_scale=3.0,
        fixed_body_axis_angle_scale=1.25,
        fixed_first_second_spring_scale=1.25,
        attach_frame_position_scales=[1.0, 2.0],
        attach_frame_tangent_scales=[1.0, 1.5],
    )

    conditions = script.build_conditions(args)

    assert [condition.condition_id for condition in conditions] == [
        "af3_axis1p25_fs1p25_fp1_ft1",
        "af3_axis1p25_fs1p25_fp1_ft1p5",
        "af3_axis1p25_fs1p25_fp2_ft1",
        "af3_axis1p25_fs1p25_fp2_ft1p5",
    ]
    assert all(condition.mode == "attach-frame-grid" for condition in conditions)
    assert conditions[-1].scales == {
        "local_attach_first_spring_scale": 3.0,
        "local_attach_first_body_axis_angle_scale": 1.25,
        "local_first_second_spring_scale": 1.25,
        "local_attach_frame_position_scale": 2.0,
        "local_attach_frame_tangent_scale": 1.5,
    }


def test_phase2_82_torque_profile_grid_conditions() -> None:
    args = SimpleNamespace(
        mode="torque-profile-grid",
        force_distributions=[
            "root_torque_segment_couples",
            "root_torque_axis_projection",
        ],
        fixed_attach_first_spring_scale=1.0,
        fixed_body_axis_angle_scale=1.0,
        fixed_first_second_spring_scale=1.0,
        fixed_attach_frame_position_scale=3.0,
        fixed_attach_frame_tangent_scale=1.5,
        torque_distribution_profiles=[
            "diffusive",
            "diffusive_sqrt",
            "diffusive_floor_0p2",
            "diffusive_floor_0p4",
            "uniform",
        ],
    )

    conditions = script.build_conditions(args)

    assert [condition.condition_id for condition in conditions[:5]] == [
        "segment_couples_diffusive_fp3_ft1p5",
        "segment_couples_diffusive_sqrt_fp3_ft1p5",
        "segment_couples_diffusive_floor_0p2_fp3_ft1p5",
        "segment_couples_diffusive_floor_0p4_fp3_ft1p5",
        "segment_couples_uniform_fp3_ft1p5",
    ]
    assert len(conditions) == 10
    assert conditions[-1].scales == {
        "force_distribution": "root_torque_axis_projection",
        "local_attach_first_spring_scale": 1.0,
        "local_attach_first_body_axis_angle_scale": 1.0,
        "local_first_second_spring_scale": 1.0,
        "local_attach_frame_position_scale": 3.0,
        "local_attach_frame_tangent_scale": 1.5,
        "torque_distribution_profile": "uniform",
    }


def test_phase2_82_axis_center_phase_summary() -> None:
    rows = [
        {"flag_id": "0", "axis_center_spin_phase_deg": "170"},
        {"flag_id": "1", "axis_center_spin_phase_deg": "10"},
        {"flag_id": "0", "axis_center_spin_phase_deg": "-170"},
        {"flag_id": "1", "axis_center_spin_phase_deg": "40"},
        {"flag_id": "0", "axis_center_spin_phase_deg": "-150"},
        {"flag_id": "1", "axis_center_spin_phase_deg": "70"},
    ]

    summary = script._axis_center_phase_summary(rows)

    assert summary["axis_center_net_abs_revolutions_min"] == pytest.approx(1.0 / 9.0)
    assert summary["axis_center_net_abs_revolutions_max"] == pytest.approx(1.0 / 6.0)
    assert summary["axis_center_net_abs_revolutions_mean"] == pytest.approx(5.0 / 36.0)
    assert summary["axis_center_direction_consistency_mean"] == 1.0
    assert summary["axis_center_direction_consistency_min"] == 1.0


def test_phase2_82_summary_row_records_fail_and_max_hook_events(
    tmp_path: Path,
) -> None:
    cfg = SimpleNamespace(
        time=SimpleNamespace(duration_s=0.5),
        dt_star=1.0e-4,
        motor_torque_Nm=2.5e-20,
        flagella=SimpleNamespace(n_flagella=3, n_beads_per_flagellum=11),
        body=SimpleNamespace(prism=SimpleNamespace(n_prism=3)),
        motor=SimpleNamespace(
            force_distribution="root_torque_segment_couples",
            local_attach_first_spring_scale=3.0,
            local_attach_first_body_axis_angle_scale=1.25,
            local_first_second_spring_scale=1.25,
            local_attach_frame_position_scale=2.0,
            local_attach_frame_tangent_scale=1.5,
            torque_distribution_profile="diffusive",
        ),
        compute_body_n_layers=lambda: 5,
    )
    condition = script.Condition(
        condition_id="af3_axis1p25_fs1p25",
        mode="first-second-grid",
        description="test",
        scales={},
    )
    rows = [
        {
            "t_s": "0.0",
            "shape_pass_nonbody": "True",
            "hook_len_rel_err_max": "0.1",
            "hook_len_rel_err_max_flag_id": "0",
            "flag_bond_rel_err_max": "0.1",
            "flag_bond_rel_err_max_flag_id": "0",
            "flag_bond_rel_err_max_bead_i": "9",
            "flag_bond_rel_err_max_bead_j": "10",
            "flag_bond_rel_err_max_len_over_b": "0.63",
            "flag_bond_rel_err_local_0_1_per_flag": "0:0.1|1:0.2|2:0.3",
            "flag_bond_rel_err_local_1_2_per_flag": "0:0.4|1:0.5|2:0.6",
            "flag_bond_rel_err_local_2_3_per_flag": "0:0.7|1:0.8|2:0.9",
            "flag_bond_rel_err_local_3_4_per_flag": "0:1.0|1:1.1|2:1.2",
            "flag_bond_rel_err_local_4_5_per_flag": "0:1.3|1:1.4|2:1.5",
        },
        {
            "t_s": "0.2",
            "shape_pass_nonbody": "False",
            "first_fail_category_nonbody": "flag",
            "hook_len_rel_err_max": "1.01",
            "hook_len_rel_err_max_flag_id": "1",
            "flag_bond_rel_err_max": "1.2",
            "flag_bond_rel_err_max_flag_id": "1",
            "flag_bond_rel_err_max_bead_i": "28",
            "flag_bond_rel_err_max_bead_j": "29",
            "flag_bond_rel_err_max_len_over_b": "1.276",
            "flag_bond_rel_err_local_0_1_per_flag": "0:0.11|1:0.21|2:0.31",
            "flag_bond_rel_err_local_1_2_per_flag": "0:0.41|1:0.51|2:0.61",
            "flag_bond_rel_err_local_2_3_per_flag": "0:0.71|1:0.81|2:0.91",
            "flag_bond_rel_err_local_3_4_per_flag": "0:1.01|1:1.11|2:1.21",
            "flag_bond_rel_err_local_4_5_per_flag": "0:1.31|1:1.41|2:1.51",
            "local_attach_first_rel_err": "0.9",
            "local_first_second_rel_err": "0.2",
            "local_attach_frame_position_rel_err": "0.3",
            "local_attach_frame_position_angle_err_deg": "4.0",
            "local_attach_frame_tangent_angle_err_deg": "5.0",
        },
        {
            "t_s": "0.5",
            "shape_pass_nonbody": "False",
            "hook_len_rel_err_max": "1.4",
            "hook_len_rel_err_max_flag_id": "2",
            "hook_len_rel_err_max_attach_body_bead_index": "8",
            "hook_len_rel_err_max_flag_first_bead_index": "37",
            "hook_len_rel_err_max_len_over_b": "0.6",
            "flag_bond_rel_err_max": "2.05",
            "flag_bond_rel_err_max_flag_id": "2",
            "flag_bond_rel_err_max_bead_i": "40",
            "flag_bond_rel_err_max_bead_j": "41",
            "flag_bond_rel_err_max_len_over_b": "1.769",
            "flag_bond_rel_err_local_0_1_per_flag": "0:0.12|1:0.22|2:0.32",
            "flag_bond_rel_err_local_1_2_per_flag": "0:0.42|1:0.52|2:0.62",
            "flag_bond_rel_err_local_2_3_per_flag": "0:0.72|1:0.82|2:0.92",
            "flag_bond_rel_err_local_3_4_per_flag": "0:1.02|1:1.12|2:1.22",
            "flag_bond_rel_err_local_4_5_per_flag": "0:1.32|1:1.42|2:1.52",
            "local_attach_frame_position_rel_err": "0.7",
            "local_attach_frame_position_angle_err_deg": "8.0",
            "local_attach_frame_tangent_angle_err_deg": "9.0",
        },
    ]

    row = script._summary_row(
        cfg,
        condition,
        tmp_path,
        rows[-1],
        rows,
        helix_summary={},
    )

    assert row["first_fail_t_s"] == "0.2"
    assert row["force_distribution"] == "root_torque_segment_couples"
    assert row["torque_distribution_profile"] == "diffusive"
    assert row["first_fail_category_nonbody"] == "flag"
    assert row["first_fail_hook_len_rel_err_max"] == "1.01"
    assert row["first_fail_hook_len_rel_err_max_flag_id"] == "1"
    assert row["first_fail_local_attach_frame_position_rel_err"] == "0.3"
    assert row["first_fail_local_attach_frame_position_angle_err_deg"] == "4.0"
    assert row["first_fail_local_attach_frame_tangent_angle_err_deg"] == "5.0"
    assert row["first_fail_flag_bond_rel_err_max"] == "1.2"
    assert row["first_fail_flag_bond_rel_err_max_flag_id"] == "1"
    assert row["first_fail_flag_bond_rel_err_max_bead_i"] == "28"
    assert row["first_fail_flag_bond_rel_err_max_bead_j"] == "29"
    assert row["first_fail_flag_bond_rel_err_max_local_bead_i"] == "2"
    assert row["first_fail_flag_bond_rel_err_max_local_bead_j"] == "3"
    assert row["first_fail_flag_bond_rel_err_max_len_over_b"] == "1.276"
    assert (
        row["first_fail_flag_bond_rel_err_local_3_4_per_flag"] == "0:1.01|1:1.11|2:1.21"
    )
    assert row["max_hook_len_rel_err_t_s"] == "0.5"
    assert row["max_hook_len_rel_err"] == "1.4"
    assert row["max_hook_len_rel_err_flag_id"] == "2"
    assert row["max_hook_len_rel_err_attach_body_bead_index"] == "8"
    assert row["max_hook_len_rel_err_flag_first_bead_index"] == "37"
    assert row["max_hook_local_attach_frame_position_rel_err"] == "0.7"
    assert row["max_hook_local_attach_frame_position_angle_err_deg"] == "8.0"
    assert row["max_hook_local_attach_frame_tangent_angle_err_deg"] == "9.0"
    assert row["max_flag_bond_rel_err_t_s"] == "0.5"
    assert row["max_flag_bond_rel_err"] == "2.05"
    assert row["max_flag_bond_rel_err_flag_id"] == "2"
    assert row["max_flag_bond_rel_err_bead_i"] == "40"
    assert row["max_flag_bond_rel_err_bead_j"] == "41"
    assert row["max_flag_bond_rel_err_local_bead_i"] == "3"
    assert row["max_flag_bond_rel_err_local_bead_j"] == "4"
    assert row["max_flag_bond_rel_err_len_over_b"] == "1.769"
    assert row["flag_bond_rel_err_max_local_bead_i"] == "3"
    assert row["flag_bond_rel_err_max_local_bead_j"] == "4"
    assert row["flag_bond_rel_err_local_3_4_per_flag"] == "0:1.02|1:1.12|2:1.22"
    assert row["max_flag_bond_rel_err_local_3_4_per_flag"] == "0:1.02|1:1.12|2:1.22"
