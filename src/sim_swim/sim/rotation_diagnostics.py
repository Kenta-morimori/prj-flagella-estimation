"""Opt-in material-motion diagnostics, integrated at every internal step.

These are best-fit axial spins of deformable labelled chains, not rotations of
PCA axes. They are diagnostics, not an additional physical acceptance gate.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np


def rigid_fit(before: np.ndarray, after: np.ndarray) -> tuple[np.ndarray, float]:
    """Return proper column-vector rotation from before to after and RMS error."""
    a = before - before.mean(0)
    b = after - after.mean(0)
    u, _, vt = np.linalg.svd(a.T @ b)
    correction = np.eye(3)
    correction[-1, -1] = np.linalg.det(vt.T @ u.T)
    rotation = vt.T @ correction @ u.T
    residual = float(np.sqrt(np.mean(np.sum((a @ rotation.T - b) ** 2, axis=1))))
    return rotation, residual


def axial_spin_increment(before: np.ndarray, after: np.ndarray) -> tuple[float, float]:
    """Finite-step best-fit rotation projected on the previous oriented PCA axis."""
    rotation, residual = rigid_fit(before, after)
    cos_angle = np.clip((np.trace(rotation) - 1) / 2, -1.0, 1.0)
    skew = (
        np.array(
            [
                rotation[2, 1] - rotation[1, 2],
                rotation[0, 2] - rotation[2, 0],
                rotation[1, 0] - rotation[0, 1],
            ]
        )
        / 2
    )
    sin_angle = np.linalg.norm(skew)
    angle = float(np.arctan2(sin_angle, cos_angle))
    if angle >= 0.25:
        return float("nan"), residual  # internal-step rotation is under-resolved
    _, singular, vt = np.linalg.svd(before - before.mean(0))
    if singular[1] < singular[0] * 1e-8:
        return float("nan"), residual  # axial phase of a straight chain is undefined
    axis = vt[0]
    if np.dot(axis, before[-1] - before[0]) < 0:
        axis = -axis
    rotvec = skew * (angle / sin_angle if sin_angle > 1e-12 else 1.0)
    return float(np.dot(rotvec, axis)), residual


class RotationDiagnosticsRecorder:
    """Integrate each step; emit compact per-flag phase/work at archive cadence."""

    def __init__(self, model, output_dir: Path, interval_s: float):
        self.model = model
        self.interval = interval_s
        self.next_t = interval_s
        self.phase = np.zeros(len(model.flagella_indices))
        self.fit_max = np.zeros_like(self.phase)
        self.valid = np.ones_like(self.phase, dtype=bool)
        self.motor_work = 0.0
        self.total_work = 0.0
        self.last_written = -1.0
        self.path = output_dir / "rotation_diagnostics.csv"
        self.file = self.path.open("w", newline="")
        fields = ["t_s", "motor_work_J", "total_force_work_J"]
        for i in range(len(self.phase)):
            fields += [f"flag{i}_phase_rad", f"flag{i}_valid", f"flag{i}_fit_rms_max_m"]
        self.writer = csv.DictWriter(self.file, fieldnames=fields)
        self.writer.writeheader()
        self.write(0.0)
        (output_dir / "rotation_diagnostics_manifest.json").write_text(
            json.dumps(
                {
                    "method": "internal-step labelled-chain rigid-fit axial spin in co-moving body frame",
                    "axis": "previous oriented PCA axis; hook first bead excluded",
                    "phase": "accumulated signed spin; never modulo 2pi",
                    "validity": "false after degenerate fit or >=0.25rad internal-step rotation; no phase reset",
                    "fit_rms": "maximum non-rigid fit residual since previous emitted sample",
                    "work": "sum F(prestep) dot displacement for all beads, in lab frame; no Brownian work separation",
                    "sampling_interval_s": interval_s,
                    "qc_gate": False,
                },
                indent=2,
            )
            + "\n"
        )

    def record(self, diag, t_s: float, *, final: bool = False):
        before, after = diag.positions_before_m, diag.positions_after_m
        body = self.model.body_indices
        rotation, _ = rigid_fit(before[body], after[body])
        aligned = (after - after[body].mean(0)) @ rotation + before[body].mean(0)
        for i, ids in enumerate(self.model.flagella_indices):
            spin, residual = axial_spin_increment(before[ids[1:]], aligned[ids[1:]])
            self.valid[i] &= np.isfinite(spin)
            self.phase[i] += spin
            self.fit_max[i] = max(self.fit_max[i], residual)
        delta = after - before
        self.motor_work += float(np.sum(diag.motor_forces * delta))
        self.total_work += float(np.sum(diag.total_forces * delta))
        if t_s + 1e-12 >= self.next_t or final:
            self.write(t_s)
            while self.next_t <= t_s + 1e-12:
                self.next_t += self.interval

    def write(self, t_s):
        if t_s <= self.last_written:
            return
        row = {
            "t_s": t_s,
            "motor_work_J": self.motor_work,
            "total_force_work_J": self.total_work,
        }
        for i in range(len(self.phase)):
            row.update(
                {
                    f"flag{i}_phase_rad": self.phase[i],
                    f"flag{i}_valid": int(self.valid[i]),
                    f"flag{i}_fit_rms_max_m": self.fit_max[i],
                }
            )
        self.writer.writerow(row)
        self.file.flush()
        self.fit_max[:] = 0
        self.last_written = t_s

    def close(self):
        self.file.close()
