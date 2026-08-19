from __future__ import annotations

from sim_swim.analysis.task_d_2015_tau_linked import _safety


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
