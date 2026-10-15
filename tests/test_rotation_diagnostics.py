import csv
from types import SimpleNamespace

import numpy as np
import pytest

from sim_swim.sim.rotation_diagnostics import (
    RotationDiagnosticsRecorder,
    axial_spin_increment,
)


def helix():
    t = np.linspace(0, 4 * np.pi, 20)
    return np.column_stack([t * 0.2, np.cos(t), np.sin(t)]) * 1e-6


def rx(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[1.0, 0, 0], [0, c, -s], [0, s, c]])


def test_axial_spin_distinguishes_rotation_and_translation():
    p = helix()
    spin, residual = axial_spin_increment(
        p, p @ rx(0.02).T + np.array([3, 2, -1]) * 1e-6
    )
    # Use the fitted axis itself for an exact axial-spin reference.
    _, _, vt = np.linalg.svd(p - p.mean(0))
    axis = vt[0]
    axis *= np.sign(np.dot(axis, p[-1] - p[0]))
    assert spin == pytest.approx(0.02 * axis[0], abs=1e-12)
    assert residual < 1e-20
    spin, _ = axial_spin_increment(p, p + 2e-6)
    assert abs(spin) < 1e-12
    invalid, _ = axial_spin_increment(p, p @ rx(0.5).T)
    assert np.isnan(invalid)


def test_recorder_removes_body_rotation_and_accumulates_before_downsample(tmp_path):
    body = np.array([[0.0, 0.0, 0.0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]) * 1e-6
    filament = helix()
    p = np.concatenate([body, filament])
    model = SimpleNamespace(
        body_indices=np.arange(4), flagella_indices=[np.arange(4, len(p))]
    )
    recorder = RotationDiagnosticsRecorder(model, tmp_path, 0.1)
    # 20 radians between output rows: cannot be recovered by archive unwrapping.
    for i in range(2000):
        angle = i * 0.01
        before = np.concatenate([body, filament @ rx(angle).T]) @ rx(i * 0.003).T
        after = (
            np.concatenate([body, filament @ rx(angle + 0.01).T])
            @ rx((i + 1) * 0.003).T
        )
        diag = SimpleNamespace(
            positions_before_m=before,
            positions_after_m=after,
            motor_forces=np.zeros_like(p),
            total_forces=np.zeros_like(p),
        )
        recorder.record(diag, (i + 1) * 0.00005, final=i == 1999)
    recorder.close()
    rows = list(csv.DictReader((tmp_path / "rotation_diagnostics.csv").open()))
    assert len(rows) == 2
    q = filament[1:]
    _, _, vt = np.linalg.svd(q - q.mean(0))
    axis = vt[0]
    axis *= np.sign(np.dot(axis, q[-1] - q[0]))
    assert float(rows[-1]["flag0_phase_rad"]) == pytest.approx(20 * axis[0], abs=1e-9)
    assert rows[-1]["flag0_valid"] == "1"
    assert float(rows[-1]["flag0_fit_rms_max_m"]) < 1e-19


def test_straight_chain_does_not_report_axial_phase():
    p = np.column_stack([np.arange(8), np.zeros(8), np.zeros(8)]) * 1e-6
    phase, _ = axial_spin_increment(p, p)
    assert np.isnan(phase)
