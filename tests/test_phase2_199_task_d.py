from __future__ import annotations

import csv

from sim_swim.analysis.task_d_2015_tau_linked import _close, _grid_key, _safety
from sim_swim.analysis.task_d_2015_tau_linked import _task_a_style_features


def test_task_d_grid_key_normalizes_floating_point_dt_representation() -> None:
    assert _grid_key(0.00010000000000000002, 1.0e-21) == _grid_key(1.0e-4, 1.0e-21)


def test_task_d_torque_comparison_distinguishes_physical_grid_values() -> None:
    assert _close(1.0e-19, 1.0e-19)
    assert not _close(1.0e-19, 2.5e-20)


def test_task_d_streams_task_a_style_axis_and_bundle_features(tmp_path) -> None:
    path = tmp_path / "step_summary.csv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "flag_helix_axis_mean_deviation_deg_max",
                "flag_helix_axis_rearward_projection_min",
                "flag_helix_bundle_radius_max_um",
                "bundle_participation_ratio",
            ],
        )
        writer.writeheader()
        writer.writerows(
            [
                {
                    "flag_helix_axis_mean_deviation_deg_max": "2.0",
                    "flag_helix_axis_rearward_projection_min": "0.99",
                    "flag_helix_bundle_radius_max_um": "1.0",
                    "bundle_participation_ratio": "0.3333333333",
                },
                {
                    "flag_helix_axis_mean_deviation_deg_max": "3.0",
                    "flag_helix_axis_rearward_projection_min": "0.95",
                    "flag_helix_bundle_radius_max_um": "1.5",
                    "bundle_participation_ratio": "0.6666666667",
                },
            ]
        )

    assert _task_a_style_features(tmp_path) == {
        "max_axis_mean_deviation_deg": 3.0,
        "min_axis_rearward_projection": 0.95,
        "max_bundle_radius_um": 1.5,
        "final_bundle_participation_ratio": 0.6666666667,
    }


def test_task_d_safety_requires_every_locked_metric_to_be_finite_and_within_limit() -> (
    None
):
    thresholds = {"body_spring_max_stretch_ratio": 0.1, "max_flag_bond_rel_err": 0.08}
    passing = {
        "status": "completed",
        "completion_pass": "True",
        "finite_pass_all": "True",
        "body_spring_max_stretch_ratio": "0.05",
        "max_flag_bond_rel_err": "0.02",
    }
    assert _safety(passing, thresholds) == (True, [])

    failing = {**passing, "max_flag_bond_rel_err": "nan"}
    passed, failures = _safety(failing, thresholds)
    assert not passed
    assert failures == ["max_flag_bond_rel_err=nonfinite"]
