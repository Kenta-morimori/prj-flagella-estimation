from __future__ import annotations

import csv
from pathlib import Path

from sim_swim.analysis.sweeps import bundling_alignment as phase27


def test_phase27_classifies_shape_fail_as_collapse() -> None:
    rows = [
        {
            "t_s": "0.4",
            "shape_pass_nonbody": "False",
            "first_fail_category_nonbody": "flag",
            "flag_helix_axis_pair_angle_deg_max": "5.0",
            "flag_helix_axis_mean_deviation_deg_max": "3.0",
        }
    ]

    assert phase27.classify_phase27_condition(rows, n_flagella=3) == "collapse"


def test_phase27_classifies_hook_wrapped_axis_aligned() -> None:
    rows = [
        {
            "t_s": "0.1",
            "shape_pass_nonbody": "True",
            "first_fail_category_nonbody": "none",
            "flag_helix_axis_pair_angle_deg_max": "5.0",
            "flag_helix_axis_mean_deviation_deg_max": "3.0",
        },
        {
            "t_s": "0.5",
            "shape_pass_nonbody": "False",
            "first_fail_category_nonbody": "hook",
            "flag_helix_axis_pair_angle_deg_max": "6.0",
            "flag_helix_axis_mean_deviation_deg_max": "4.0",
        },
    ]

    assert (
        phase27.classify_phase27_condition(rows, n_flagella=3)
        == "hook_wrapped_axis_aligned"
    )


def test_phase27_classifies_axis_aligned_stable() -> None:
    rows = [
        {
            "t_s": "0.1",
            "shape_pass_nonbody": "True",
            "first_fail_category_nonbody": "none",
            "flag_helix_axis_pair_angle_deg_max": "5.0",
            "flag_helix_axis_mean_deviation_deg_max": "3.0",
        },
        {
            "t_s": "0.5",
            "shape_pass_nonbody": "True",
            "first_fail_category_nonbody": "none",
            "flag_helix_axis_pair_angle_deg_max": "10.0",
            "flag_helix_axis_mean_deviation_deg_max": "5.0",
        },
    ]

    assert (
        phase27.classify_phase27_condition(rows, n_flagella=3) == "axis_aligned_stable"
    )


def test_phase27_classifies_one_outlier_as_axis_not_aligned() -> None:
    rows = [
        {
            "t_s": "0.1",
            "shape_pass_nonbody": "True",
            "first_fail_category_nonbody": "none",
            "flag_helix_axis_pair_angle_deg_max": "5.0",
            "flag_helix_axis_mean_deviation_deg_max": "3.0",
        },
        {
            "t_s": "0.5",
            "shape_pass_nonbody": "True",
            "first_fail_category_nonbody": "none",
            "flag_helix_axis_pair_angle_deg_max": "45.0",
            "flag_helix_axis_mean_deviation_deg_max": "30.0",
        },
    ]

    assert phase27.classify_phase27_condition(rows, n_flagella=3) == "axis_not_aligned"


def test_phase27_temporal_stability_ignores_initial_20_percent() -> None:
    rows = [
        {
            "t_s": "0.0",
            "shape_pass_nonbody": "True",
            "first_fail_category_nonbody": "none",
            "flag_helix_axis_pair_angle_deg_max": "60.0",
            "flag_helix_axis_mean_deviation_deg_max": "40.0",
        },
        {
            "t_s": "0.2",
            "shape_pass_nonbody": "True",
            "first_fail_category_nonbody": "none",
            "flag_helix_axis_pair_angle_deg_max": "10.0",
            "flag_helix_axis_mean_deviation_deg_max": "5.0",
        },
        {
            "t_s": "0.5",
            "shape_pass_nonbody": "True",
            "first_fail_category_nonbody": "none",
            "flag_helix_axis_pair_angle_deg_max": "8.0",
            "flag_helix_axis_mean_deviation_deg_max": "4.0",
        },
    ]

    assert (
        phase27.classify_phase27_condition(rows, n_flagella=3) == "axis_aligned_stable"
    )


def test_phase27_net_abs_helix_spin_revolutions_uses_phase_delta() -> None:
    rows = [
        {"flag_helix_spin_phase_deg": "10.0"},
        {"flag_helix_spin_phase_deg": "370.0"},
    ]

    assert phase27._net_abs_helix_spin_revolutions(rows) == 1.0


def test_phase27_build_config_sets_initial_helix_axis_angle() -> None:
    raw_cfg = phase27._load_yaml(
        Path(__file__).resolve().parents[1] / "conf/sim_swim.yaml"
    )

    cfg = phase27._build_config(
        raw_cfg,
        initial_helix_axis_from_rear_deg=0.0,
        n_flagella=3,
        torque_Nm=0.5e-20,
        duration_s=0.02,
        dt_star=1.0e-4,
        overrides=[],
    )

    assert cfg.flagella.initial_helix_axis_from_rear_deg == 0.0
    assert cfg.flagella.n_flagella == 3
    assert cfg.motor.torque_Nm == 0.5e-20


def test_phase27_plot_flag_helix_axis_timeseries_outputs_png(tmp_path: Path) -> None:
    axis_csv = tmp_path / "flag_helix_axis_diagnostics.csv"
    with axis_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "step",
                "t_s",
                "flag_id",
                "flag_helix_axis_vs_rear_angle_deg",
                "flag_helix_axis_rearward_projection",
                "flag_helix_axis_fit_r2",
                "axis_origin_x_um",
                "axis_origin_y_um",
                "axis_origin_z_um",
                "axis_dir_x",
                "axis_dir_y",
                "axis_dir_z",
            ],
        )
        writer.writeheader()
        for step, t_s in enumerate([0.0, 0.1, 0.2]):
            rows = [
                (0, 5.0, 1.0, 0.0, 0.0),
                (1, 6.0, 0.99, 0.05, 0.0),
                (2, 80.0, 0.2, 0.98, 0.0),
            ]
            for flag_id, rear_angle, axis_x, axis_y, axis_z in rows:
                writer.writerow(
                    {
                        "step": step,
                        "t_s": t_s,
                        "flag_id": flag_id,
                        "flag_helix_axis_vs_rear_angle_deg": rear_angle,
                        "flag_helix_axis_rearward_projection": 1.0,
                        "flag_helix_axis_fit_r2": 0.9,
                        "axis_origin_x_um": 0.0,
                        "axis_origin_y_um": 0.0,
                        "axis_origin_z_um": 0.0,
                        "axis_dir_x": axis_x,
                        "axis_dir_y": axis_y,
                        "axis_dir_z": axis_z,
                    }
                )

    output_path = tmp_path / "flag_helix_axis_angles_timeseries.png"

    phase27.plot_flag_helix_axis_timeseries(
        axis_csv_path=axis_csv,
        output_path=output_path,
    )

    assert output_path.is_file()
    assert output_path.stat().st_size > 0
