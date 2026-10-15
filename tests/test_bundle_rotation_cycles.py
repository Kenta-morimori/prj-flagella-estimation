import numpy as np
import pytest

from sim_swim.analysis.bundle_rotation import (
    completed_rotation_cycles,
    observation_windows,
)


def rows(direction=1):
    return [
        {
            "condition_id": "synthetic",
            "t_s": float(t),
            "flag0_phase_rad": direction * 4 * np.pi * t,
            "flag0_valid": True,
            "bundle_radius_mean_um": 1.0 - 0.1 * t,
            "helix_pair_angle_mean_deg": 20.0 - t,
            "bundle_rear_angle_deg": 2.0,
        }
        for t in np.linspace(0.0, 1.25, 126)
    ]


@pytest.mark.parametrize("direction", [1, -1])
def test_full_cycles_preserve_sign_and_measure_boundary_drift(direction):
    cycles, status = completed_rotation_cycles(rows(direction), 1)
    assert status.startswith("available") and len(cycles) == 2
    for cycle in cycles:
        assert cycle["period_s"] == pytest.approx(0.5)
        assert cycle["flag0_net_turns"] == pytest.approx(direction)
        assert cycle["bundle_radius_mean_um_boundary_change"] == pytest.approx(-0.05)
    assert cycles[-1]["end_s"] == pytest.approx(
        1.0
    )  # remaining quarter second is not a cycle


def test_missing_phase_and_reversal_are_not_completed_cycles():
    a = rows()
    a[20]["flag0_valid"] = False
    assert completed_rotation_cycles(a, 1)[0] == []
    a = rows()
    a[40]["flag0_phase_rad"] = a[39]["flag0_phase_rad"] - 0.1
    assert completed_rotation_cycles(a, 1)[0] == []
    a = rows()
    a[10]["flag0_phase_rad"] = float("nan")
    assert completed_rotation_cycles(a, 1)[0] == []


def test_longer_run_has_actual_terminal_windows_and_no_duplicates():
    windows = observation_windows(1.999999999998)
    assert (1.5, 2.0) in windows and (1.75, 2.0) in windows
    assert len(set(windows)) == len(windows)
    assert observation_windows(0.04) == [(0.0, 0.04)]


def test_accumulated_turns_do_not_imply_resolved_cycle_geometry():
    a = rows()[::50]
    assert completed_rotation_cycles(a, 1)[0] == []
    assert "sampling" in completed_rotation_cycles(a, 1)[1]
