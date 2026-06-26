"""Hook attach local-frame helpers."""

from __future__ import annotations

import numpy as np

from sim_swim.model.types import SimModel


def _unit_vector(vec: np.ndarray) -> np.ndarray:
    norm = float(np.linalg.norm(vec))
    if norm <= 1e-18:
        return np.zeros(3, dtype=float)
    return vec / norm


def _fallback_perpendicular(axis: np.ndarray) -> np.ndarray:
    ref = np.array([1.0, 0.0, 0.0], dtype=float)
    if abs(float(np.dot(_unit_vector(axis), ref))) > 0.9:
        ref = np.array([0.0, 1.0, 0.0], dtype=float)
    out = _unit_vector(np.cross(axis, ref))
    if float(np.linalg.norm(out)) <= 1e-18:
        return np.array([0.0, 0.0, 1.0], dtype=float)
    return out


def hook_attach_layer_indices(model: SimModel) -> np.ndarray:
    """Return the body layer index for each hook attach bead."""

    rows: list[int] = []
    for attach_raw, _first_raw, _second_raw in model.hook_triplets.astype(
        int, copy=False
    ):
        attach = int(attach_raw)
        layer_idx = -1
        for idx, layer in enumerate(model.body_layer_indices):
            if np.any(layer.astype(int, copy=False) == attach):
                layer_idx = int(idx)
                break
        rows.append(layer_idx)
    return np.asarray(rows, dtype=int)


def hook_attach_frame(
    positions_m: np.ndarray,
    model: SimModel,
    attach_idx: int,
    layer_idx: int,
    body_axis_unit: np.ndarray,
) -> np.ndarray:
    """Build an orthonormal frame tied to one body surface attach point."""

    axis = _unit_vector(np.asarray(body_axis_unit, dtype=float))
    if float(np.linalg.norm(axis)) <= 1e-18:
        axis = np.array([1.0, 0.0, 0.0], dtype=float)

    radial = np.zeros(3, dtype=float)
    if 0 <= int(layer_idx) < len(model.body_layer_indices):
        layer = model.body_layer_indices[int(layer_idx)].astype(int, copy=False)
        if layer.size > 0:
            centroid = np.mean(positions_m[layer], axis=0)
            radial = positions_m[int(attach_idx)] - centroid
    radial = radial - float(np.dot(radial, axis)) * axis
    radial = _unit_vector(radial)
    if float(np.linalg.norm(radial)) <= 1e-18:
        radial = _fallback_perpendicular(axis)

    tangent = _unit_vector(np.cross(axis, radial))
    if float(np.linalg.norm(tangent)) <= 1e-18:
        radial = _fallback_perpendicular(axis)
        tangent = _unit_vector(np.cross(axis, radial))
    return np.column_stack([axis, radial, tangent])


def hook_frame_local_vectors(
    positions_m: np.ndarray,
    model: SimModel,
    attach_layer_indices: np.ndarray,
    body_axis_unit: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Record hook-root vectors in each attach point's local frame."""

    n_hook = int(model.hook_triplets.shape[0])
    attach_first = np.zeros((n_hook, 3), dtype=float)
    first_second = np.zeros((n_hook, 3), dtype=float)
    for row, (attach_raw, first_raw, second_raw) in enumerate(
        model.hook_triplets.astype(int, copy=False)
    ):
        attach = int(attach_raw)
        first = int(first_raw)
        second = int(second_raw)
        layer_idx = (
            int(attach_layer_indices[row])
            if row < int(attach_layer_indices.shape[0])
            else -1
        )
        frame = hook_attach_frame(
            positions_m=positions_m,
            model=model,
            attach_idx=attach,
            layer_idx=layer_idx,
            body_axis_unit=body_axis_unit,
        )
        attach_first[row] = frame.T @ (positions_m[first] - positions_m[attach])
        first_second[row] = frame.T @ (positions_m[second] - positions_m[first])
    return attach_first, first_second


def hook_frame_target_vectors(
    positions_m: np.ndarray,
    model: SimModel,
    attach_layer_indices: np.ndarray,
    attach_first_local_m: np.ndarray,
    first_second_local_m: np.ndarray,
    body_axis_unit: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Transform stored local hook vectors into current world coordinates."""

    n_hook = int(model.hook_triplets.shape[0])
    attach_first = np.zeros((n_hook, 3), dtype=float)
    first_second = np.zeros((n_hook, 3), dtype=float)
    for row, (attach_raw, _first_raw, _second_raw) in enumerate(
        model.hook_triplets.astype(int, copy=False)
    ):
        attach = int(attach_raw)
        layer_idx = (
            int(attach_layer_indices[row])
            if row < int(attach_layer_indices.shape[0])
            else -1
        )
        frame = hook_attach_frame(
            positions_m=positions_m,
            model=model,
            attach_idx=attach,
            layer_idx=layer_idx,
            body_axis_unit=body_axis_unit,
        )
        if row < int(attach_first_local_m.shape[0]):
            attach_first[row] = frame @ attach_first_local_m[row]
        if row < int(first_second_local_m.shape[0]):
            first_second[row] = frame @ first_second_local_m[row]
    return attach_first, first_second
