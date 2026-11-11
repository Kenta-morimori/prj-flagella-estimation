"""Opt-in effective motor: analytic virtual work of a relative frame phase.

This is not a material-twist rod. The input frame fits the labelled local
body support; the output fits the labelled chain and minimally aligns its
reference axis with the current oriented PCA axis. For their relative
quaternion, exp(i chi)=(w+i*z)^2/(w*w+z*z), and F=T*grad(chi).

If a is the raw fitted flag axis, e the PCA axis, n the input axis, then
g=(n+e)/(1+n.e), k=(a+e)/(1+a.e), h=(a cross e)/(1+a.e) give
dchi=k.dtheta_fit + (g cross e-h).de - g.dtheta_body.
Both the proper-fit and PCA terms below are analytic pullbacks of this
identity. A minimum-norm torque couple alone is NOT this gradient on a
deformable chain. Numerical differentiation is used only by tests.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.spatial.transform import Rotation

from sim_swim.dynamics.forces import MotorForceDiagnostics

_CONDITION_EPS = 1e-8  # numerical nonuniqueness/singularity guard, not a QC fit


def _unit(v: np.ndarray) -> np.ndarray:
    length = float(np.linalg.norm(v))
    if not np.isfinite(length) or length <= 0:
        raise RuntimeError("Phase drive: undefined reference axis")
    return v / length


def _skew(v: np.ndarray) -> np.ndarray:
    x, y, z = v
    return np.array([[0.0, -z, y], [z, 0.0, -x], [-y, x, 0.0]])


def _align(a: np.ndarray, e: np.ndarray) -> np.ndarray:
    denominator = 1 + float(a @ e)
    if denominator <= _CONDITION_EPS:
        raise RuntimeError("Phase drive: antiparallel frame alignment")
    cross = _skew(np.cross(a, e))
    return np.eye(3) + cross + cross @ cross / denominator


@dataclass
class _Fit:
    reference: np.ndarray
    length: float

    @classmethod
    def from_points(cls, points: np.ndarray) -> _Fit:
        if len(points) < 3:
            raise RuntimeError(
                "Phase drive: material support needs at least three points"
            )
        centered = points - points.mean(0)
        length = float(np.max(np.linalg.norm(centered, axis=1)))
        if not np.isfinite(length) or length <= 0:
            raise RuntimeError("Phase drive: degenerate reference support")
        return cls(centered / length, length)

    def evaluate(self, points: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
        current = (points - points.mean(0)) / self.length
        cross = current.T @ self.reference
        u, _, vt = np.linalg.svd(cross)
        sign = np.diag([1.0, 1.0, 1.0 if np.linalg.det(u @ vt) >= 0 else -1.0])
        rotation = u @ sign @ vt
        stretch = cross @ rotation.T
        stretch = (stretch + stretch.T) / 2
        jacobian = np.trace(stretch) * np.eye(3) - stretch
        spectrum = np.linalg.eigvalsh(jacobian)
        if spectrum[-1] <= 0 or spectrum[0] <= _CONDITION_EPS * spectrum[-1]:
            raise RuntimeError(
                f"Phase drive: nonunique proper material fit, spectrum={spectrum.tolist()}"
            )
        return rotation, jacobian, float(spectrum[-1] / spectrum[0])

    def pullback(
        self, rotation: np.ndarray, jacobian: np.ndarray, couple: np.ndarray
    ) -> np.ndarray:
        adjoint = np.linalg.solve(jacobian, couple)
        return np.cross(adjoint, self.reference @ rotation.T) / self.length


def _pca(points: np.ndarray, length: float):
    centered = (points - points.mean(0)) / length
    values, vectors = np.linalg.eigh(centered.T @ centered)
    if values[-1] <= 0:
        raise RuntimeError("Phase drive: degenerate output PCA axis")
    gap = (values[-1] - values[-2]) / values[-1]
    direction = points[-1] - points[0]
    axis = vectors[:, -1]
    projection = float(axis @ direction)
    if gap <= _CONDITION_EPS or abs(projection) <= _CONDITION_EPS * np.linalg.norm(
        direction
    ):
        raise RuntimeError(
            f"Phase drive: degenerate or unoriented output PCA axis, gap={gap:.6g}, projection={projection:.6g}"
        )
    if projection < 0:
        axis = -axis
    return axis, gap, centered, values, vectors


def _pca_pullback(centered, values, vectors, axis, covector, length):
    # de=sum_j v_j (v_j^T dC e)/(lambda_max-lambda_j).
    other = vectors[:, :2]
    adjoint = other @ ((other.T @ covector) / (values[-1] - values[:2]))
    return (
        np.outer(centered @ axis, adjoint) + np.outer(centered @ adjoint, axis)
    ) / length


@dataclass
class _Reference:
    flag_ids: np.ndarray
    body_ids: np.ndarray
    flag_fit: _Fit
    body_fit: _Fit
    input_frame: np.ndarray
    output_frame: np.ndarray


class PhaseWorkDrive:
    """Fixed material templates; no extra dynamical state or restoring energy."""

    def __init__(self, positions, flagella_indices, body_supports, input_axes):
        self.references = []
        for ids, support, normal in zip(
            flagella_indices, body_supports, input_axes, strict=True
        ):
            ids, support = np.asarray(ids, int), np.asarray(support, int)
            flag_fit = _Fit.from_points(positions[ids])
            body_fit = _Fit.from_points(positions[support])
            flag_fit.evaluate(positions[ids])
            body_fit.evaluate(positions[support])
            n = _unit(np.asarray(normal, float))
            e = _pca(positions[ids], flag_fit.length)[0]
            # A transverse gauge, not a physical preferred direction.
            x = e - (e @ n) * n
            if np.linalg.norm(x) < _CONDITION_EPS:
                x = np.eye(3)[np.argmin(abs(n))]
                x = x - (x @ n) * n
            x = _unit(x)
            body_frame = np.column_stack([x, np.cross(n, x), n])
            self.references.append(
                _Reference(
                    ids,
                    support,
                    flag_fit,
                    body_fit,
                    body_frame,
                    _align(n, e) @ body_frame,
                )
            )
        self.work_total_J = 0.0  # diagnostic integral only; never used in forces
        self._before_phases: dict[int, complex] = {}

    def _evaluate(self, positions, ref):
        rf, jf, cf = ref.flag_fit.evaluate(positions[ref.flag_ids])
        rb, jb, cb = ref.body_fit.evaluate(positions[ref.body_ids])
        e, gap, centered, values, vectors = _pca(
            positions[ref.flag_ids], ref.flag_fit.length
        )
        n = rb @ ref.input_frame[:, 2]
        a = rf @ ref.output_frame[:, 2]
        joint, align = 1 + float(n @ e), 1 + float(a @ e)
        if joint <= _CONDITION_EPS:
            raise RuntimeError(
                f"Phase drive: antiparallel input/output axes, 1+dot={joint:.6g}"
            )
        frame = _align(a, e) @ rf @ ref.output_frame
        relative = (rb @ ref.input_frame).T @ frame
        quat = Rotation.from_matrix(relative).as_quat()
        z, w = quat[2:]
        denominator = w * w + z * z
        if denominator <= _CONDITION_EPS / 2:
            raise RuntimeError("Phase drive: undefined relative spin phase")
        phase = complex((w + 1j * z) ** 2 / denominator)
        return (
            rf,
            jf,
            rb,
            jb,
            e,
            n,
            a,
            joint,
            align,
            phase,
            gap,
            centered,
            values,
            vectors,
            max(cf, cb),
        )

    def phases(self, positions) -> np.ndarray:
        """Wrapped phase values for independent tests/archive diagnostics."""
        return np.array([self._evaluate(positions, ref)[9] for ref in self.references])

    def forces(self, positions, torque_per_flag):
        if not np.isfinite(positions).all():
            raise RuntimeError("Phase drive: nonfinite positions")
        if len(torque_per_flag) != len(self.references):
            raise ValueError("Phase drive requires one signed torque per filament")
        forces = np.zeros_like(positions, dtype=float)
        self._before_phases = {}
        records = []
        flag_torques, body_torques, flag_forces, body_forces, supports = (
            [],
            [],
            [],
            [],
            [],
        )
        for i, ref in enumerate(self.references):
            torque = float(torque_per_flag[i])
            if not np.isfinite(torque):
                raise RuntimeError(f"Phase drive: nonfinite torque for flag {i}")
            if torque == 0:
                continue
            (
                rf,
                jf,
                rb,
                jb,
                e,
                n,
                a,
                joint,
                align,
                phase,
                gap,
                centered,
                values,
                vectors,
                condition,
            ) = self._evaluate(positions, ref)
            g, k = (n + e) / joint, (a + e) / align
            h = np.cross(a, e) / align
            ff = ref.flag_fit.pullback(rf, jf, torque * k)
            ff += _pca_pullback(
                centered,
                values,
                vectors,
                e,
                torque * (np.cross(g, e) - h),
                ref.flag_fit.length,
            )
            fb = ref.body_fit.pullback(rb, jb, -torque * g)
            forces[ref.flag_ids] += ff
            forces[ref.body_ids] += fb
            target = torque * g
            rflag = positions[ref.flag_ids] - positions[ref.flag_ids].mean(0)
            rbody = positions[ref.body_ids] - positions[ref.body_ids].mean(0)
            tf, tb = np.cross(rflag, ff).sum(0), np.cross(rbody, fb).sum(0)
            records.append(
                {
                    "target_vector_error_ratio_max": np.linalg.norm(tf - target)
                    / abs(torque),
                    "output_projection_error_ratio_max": abs(tf @ e - torque)
                    / abs(torque),
                    "input_projection_error_ratio_max": abs((-tb) @ n - torque)
                    / abs(torque),
                    "transverse_ratio_max": np.linalg.norm(tf - (tf @ e) * e)
                    / abs(torque),
                    "force_norm_max_N": np.linalg.norm(ff),
                    "input_output_angle_deg_max": np.degrees(
                        np.arccos(np.clip(n @ e, -1, 1))
                    ),
                    "gain_max": np.linalg.norm(g),
                    "phase_fit_condition_max": condition,
                    "phase_pca_gap_min": gap,
                    "phase_alignment_one_plus_dot_min": align,
                    "phase_joint_one_plus_dot_min": joint,
                }
            )
            self._before_phases[i] = phase
            flag_torques.append(np.linalg.norm(tf))
            body_torques.append(np.linalg.norm(tb))
            flag_forces.append(np.linalg.norm(ff, axis=1).mean())
            body_forces.append(np.linalg.norm(fb, axis=1).mean())
            supports.append(len(ref.body_ids))
        metrics = {"motor_drive_contract_version": 2.0}
        # The historical metric is ||tau-T*e||/T for pure D, not the M error.
        metrics["motor_drive_axial_error_ratio_max"] = float("nan")
        if records:
            for name in records[0]:
                reduce = min if name.endswith("_min") else max
                metrics["motor_drive_" + name] = float(reduce(r[name] for r in records))

        def mean(values):
            return float(np.mean(values)) if values else 0.0

        return forces, MotorForceDiagnostics(
            Ta_norm_mean=mean(flag_torques),
            Tb_norm_mean=mean(body_torques),
            Fa_norm_mean=mean(flag_forces),
            Fb_norm_mean=mean(body_forces),
            reaction_support_bead_counts=tuple(supports),
            reaction_fallback_used=False,
            drive_metrics=tuple(metrics.items()),
        )

    def record_displacement(self, before, after, torques, forces):
        """Work telemetry across the actual step; no feedback into dynamics."""
        work, increment = 0.0, 0.0
        if not np.array_equal(before, after):
            for i, old in self._before_phases.items():
                new = self._evaluate(after, self.references[i])[9]
                delta = float(np.angle(new / old))
                if abs(delta) >= 0.25:
                    raise RuntimeError(
                        f"Phase drive: under-resolved internal phase step, flag={i}, increment={delta:.6g} rad"
                    )
                work += float(torques[i]) * delta
                increment = max(increment, abs(delta))
        linear = float(np.sum(forces * (after - before)))
        self.work_total_J += work
        return {
            "motor_drive_phase_work_step_J": work,
            "motor_drive_phase_linear_work_step_J": linear,
            "motor_drive_phase_work_total_J": self.work_total_J,
            "motor_drive_phase_increment_abs_max_rad": increment,
            "motor_drive_phase_work_difference_abs_J": abs(work - linear),
        }
