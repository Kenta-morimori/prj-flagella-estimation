"""Opt-in axis-based root transition with fixed, local reference geometry."""

from __future__ import annotations

import math

import numpy as np


def _frame(points: np.ndarray) -> np.ndarray:
    x = points[1] - points[0]
    length = float(np.linalg.norm(x))
    z = np.cross(x, points[2] - points[0])
    area = float(np.linalg.norm(z))
    if length <= 0 or area <= 1e-10 * length**2:
        raise ValueError("Root transition requires noncollinear reference beads")
    x, z = x / length, z / area
    return np.column_stack((x, np.cross(z, x), z))


def _rotate(axis: np.ndarray, angle: float) -> np.ndarray:
    x, y, z = axis
    K = np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]])
    return np.eye(3) + math.sin(angle) * K + (1 - math.cos(angle)) * (K @ K)


def place_axis_transition(
    body_um: np.ndarray,
    chains_um: list[np.ndarray],
    attach_ids: np.ndarray,
    body_axis: np.ndarray,
    hook_length_um: float,
    b_um: float,
    phase_offsets_deg: tuple[float, ...],
    bead_diameter_um: float,
) -> dict:
    """Construct a tapered prefix; preserve every bond length and distal helix.

    Radii 0/.035/.040 and 10 degrees were fixed by geometric and passive
    tangent-stability audits before integration. They are not runtime fitting
    controls. The source helix supplies its actual screw radius/phase increment.
    """
    rear = -body_axis / np.linalg.norm(body_axis)
    records = []
    if len(phase_offsets_deg) != len(chains_um):
        raise ValueError("Root transition phase count mismatch")
    for points, attach_id, offset in zip(chains_um, attach_ids, phase_offsets_deg):
        if len(points) < 7:
            raise ValueError("Root transition requires at least seven beads")
        old = points.copy()
        step_rotation = _frame(old[1:4]) @ _frame(old[:3]).T
        eigenvalues, vectors = np.linalg.eig(step_rotation)
        axis = vectors[:, np.argmin(abs(eigenvalues - 1))].real
        axis /= np.linalg.norm(axis)
        if axis @ (old[-1] - old[0]) < 0:
            axis = -axis
        shift = old[1] - step_rotation @ old[0]
        center = np.linalg.lstsq(
            np.vstack((np.eye(3) - step_rotation, axis)),
            np.r_[shift - (shift @ axis) * axis, axis @ old[0]],
            rcond=None,
        )[0]
        residual = np.max(
            np.linalg.norm(old[1:] - (old[:-1] @ step_rotation.T + shift), axis=1)
        )
        if residual > 1e-8 * b_um:
            raise ValueError("Root transition reference is not a uniform helix")
        rho0 = old[0] - center
        rho1 = old[1] - center
        rho1 -= (rho1 @ axis) * axis
        radius = float(np.linalg.norm(rho0))
        phase_step = math.atan2(float(axis @ np.cross(rho0, rho1)), float(rho0 @ rho1))
        cross = np.cross(axis, rear)
        sin_angle = float(np.linalg.norm(cross))
        cos_angle = float(np.clip(axis @ rear, -1, 1))
        if sin_angle <= 1e-12:
            if cos_angle < 0:
                raise ValueError("Root transition antiparallel reference axis")
            align = np.eye(3)
        else:
            align = _rotate(cross / sin_angle, math.atan2(sin_angle, cos_angle))
        rho = align @ rho0
        rho -= (rho @ rear) * rear
        rho /= np.linalg.norm(rho)
        radial = body_um[int(attach_id)] - body_um.mean(axis=0)
        radial -= (radial @ rear) * rear
        radial /= np.linalg.norm(radial)
        angle = math.radians(10)
        anchor = body_um[int(attach_id)] + hook_length_um * (
            math.cos(angle) * radial + math.sin(angle) * rear
        )
        radii = np.r_[0, 0.035 * b_um, 0.040 * b_um, np.full(len(points) - 3, radius)]
        lateral = np.array(
            [
                (_rotate(rear, i * phase_step + math.radians(offset)) @ rho) * r
                for i, r in enumerate(radii)
            ]
        )
        lengths = np.linalg.norm(np.diff(old, axis=0), axis=1)
        squared = lengths**2 - np.sum(np.diff(lateral, axis=0) ** 2, axis=1)
        if np.any(squared <= 0):
            raise ValueError("Root transition cannot preserve reference bond lengths")
        axial = np.r_[0, np.cumsum(np.sqrt(squared))]
        points[:] = anchor + np.outer(axial, rear) + lateral
        records.append(
            {
                "axis_anchor_um": anchor.tolist(),
                "axis_unit": rear.tolist(),
                "radii_um": radii.tolist(),
                "phase_offset_deg": float(offset),
                "source_screw_radius_um": radius,
                "source_screw_residual_um": float(residual),
                "reference_prefix_start_indices": [0, 1, 2],
            }
        )
    # Bound the distance for every continuous spin phase. A cylinder contains
    # the body-centre hull; the one-Lipschitz distance function bounds errors
    # between samples along each line segment. No spin/time integration occurs.
    center_body = body_um.mean(axis=0)
    axial_body = (body_um - center_body) @ body_axis
    radial_body = (body_um - center_body) - np.outer(axial_body, body_axis)
    body_radius = float(np.max(np.linalg.norm(radial_body, axis=1)))
    interpolation = np.linspace(0, 1, 1001)
    lower_bounds = []
    for points, record in zip(chains_um, records):
        anchor = np.asarray(record["axis_anchor_um"])
        offset = anchor - center_body
        radius_axis = float(np.linalg.norm(offset - (offset @ body_axis) * body_axis))
        for a, b in zip(points[:-1], points[1:]):
            sample = (1 - interpolation[:, None]) * a + interpolation[:, None] * b
            axial = (sample - center_body) @ body_axis
            relative = sample - anchor
            lateral = relative - np.outer(relative @ body_axis, body_axis)
            phase_radius = np.linalg.norm(lateral, axis=1)
            distance = np.hypot(
                np.maximum(radius_axis - phase_radius - body_radius, 0),
                np.maximum(
                    np.maximum(axial_body.min() - axial, axial - axial_body.max()), 0
                ),
            )
            lower_bounds.append(float(distance.min() - np.linalg.norm(b - a) / 2000))
    bound = min(lower_bounds, default=math.inf)
    if bound < bead_diameter_um:
        raise ValueError("Root transition lacks full-phase body clearance")
    return {
        "all_phase_body_clearance_lower_bound_um": bound,
        "body_clearance_comparison_um": bead_diameter_um,
        "axis_definition": "reference_screw_axis_aligned_with_body_rear",
        "method": "axis_transition_rear_v1",
        "filaments": records,
        "hook_outward_to_rear_deg": 10.0,
        "natural_geometry_changed": True,
        "reference_scope": "proximal bend/torsion rows starting at beads 0..2 only",
        "already_dense_bundle": False,
    }
