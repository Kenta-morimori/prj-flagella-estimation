"""Conservative attachment-frame potential, including body-frame derivatives."""

from __future__ import annotations

import numpy as np

from sim_swim.model.types import SimModel


def _normalize(vector: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    norm = float(np.linalg.norm(vector))
    if not np.isfinite(norm) or norm <= 1e-18:
        raise RuntimeError("Attachment-frame energy gradient: degenerate frame")
    unit = vector / norm
    return unit, (np.eye(3) - np.outer(unit, unit)) / norm


def attachment_frame_energy_forces(
    positions_m: np.ndarray,
    model: SimModel,
    layer_indices: np.ndarray,
    local_vectors: tuple[np.ndarray, np.ndarray],
    rest_lengths: tuple[np.ndarray, np.ndarray],
    stiffness: tuple[float, float],
) -> tuple[float, np.ndarray]:
    """Return E and -grad(E); no frozen-frame approximation or world-axis fallback."""
    forces = np.zeros_like(positions_m)
    energy = 0.0
    if not any(k > 0.0 for k in stiffness):
        return energy, forces
    if len(model.body_layer_indices) < 2:
        raise RuntimeError("Attachment-frame energy gradient requires body layers")
    start, end = (np.asarray(model.body_layer_indices[i], dtype=int) for i in (0, -1))
    axis, axis_jac = _normalize(positions_m[end].mean(0) - positions_m[start].mean(0))
    for row, (attach, first, second) in enumerate(model.hook_triplets):
        layer_index = int(layer_indices[row])
        if not 0 <= layer_index < len(model.body_layer_indices):
            raise RuntimeError("Attachment-frame energy gradient requires attach layer")
        layer = np.asarray(model.body_layer_indices[layer_index], dtype=int)
        offset = positions_m[attach] - positions_m[layer].mean(0)
        radial, radial_jac = _normalize(offset - np.dot(offset, axis) * axis)
        tangent, tangent_jac = _normalize(np.cross(axis, radial))
        frame = np.column_stack((axis, radial, tangent))
        for term, (left, right) in enumerate(((attach, first), (first, second))):
            k = stiffness[term]
            if k <= 0.0:
                continue
            local = local_vectors[term][row]
            delta = positions_m[right] - positions_m[left] - frame @ local
            coefficient = k / max(float(rest_lengths[term][row]), 1e-18) ** 2
            energy += 0.5 * coefficient * float(delta @ delta)
            q = coefficient * delta
            forces[right] -= q
            forces[left] += q
            # Pull back the force covector on the target through the frame.
            gt = tangent_jac @ (local[2] * q)
            ga = local[0] * q + np.cross(radial, gt)
            gr = local[1] * q + np.cross(gt, axis)
            gv = radial_jac @ gr
            gu = gv - axis * np.dot(axis, gv)
            ga -= np.dot(offset, axis) * gv + np.dot(gv, axis) * offset
            gd = axis_jac @ ga
            forces[attach] += gu
            forces[layer] -= gu / len(layer)
            forces[end] += gd / len(end)
            forces[start] -= gd / len(start)
    return energy, forces
