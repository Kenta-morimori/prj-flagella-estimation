"""Validate and load an opt-in, archived initial bundle without changing physics."""

from __future__ import annotations

import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np


def segment_distance_squared(p, q, r, s):
    """Exact finite-segment distances in normalized model coordinates.

    Check the unconstrained interior solution and all four boundary candidates.
    Scale the parallel tolerance by segment lengths rather than physical units.
    """
    u, v, w = q - p, s - r, p - r
    a = np.sum(u * u, axis=-1)
    b = np.sum(u * v, axis=-1)
    c = np.sum(v * v, axis=-1)
    d = np.sum(u * w, axis=-1)
    e = np.sum(v * w, axis=-1)
    if np.any(a <= 0) or np.any(c <= 0):
        raise ValueError("Initial bundle contains a zero-length segment")
    denominator = a * c - b * b
    parallel = denominator <= 1e-12 * a * c
    safe = np.where(parallel, 1.0, denominator)
    first = (b * e - c * d) / safe
    second = (a * e - b * d) / safe
    interior = np.sum((w + first[:, None] * u - second[:, None] * v) ** 2, axis=-1)
    interior[parallel | (first < 0) | (first > 1) | (second < 0) | (second > 1)] = (
        np.inf
    )
    one = np.sum((w - np.clip(e / c, 0, 1)[:, None] * v) ** 2, axis=-1)
    wq = q - r
    two = np.sum(
        (wq - np.clip(np.sum(wq * v, axis=-1) / c, 0, 1)[:, None] * v) ** 2, axis=-1
    )
    three = np.sum((r - p - np.clip(-d / a, 0, 1)[:, None] * u) ** 2, axis=-1)
    ws = s - p
    four = np.sum(
        (ws - np.clip(np.sum(ws * u, axis=-1) / a, 0, 1)[:, None] * u) ** 2, axis=-1
    )
    return np.minimum.reduce([interior, one, two, three, four])


def load_initial_bundle(cfg, body_um, original_chains_um, attach_ids, body_edges):
    """Load a fixed state with verified topology, bond lengths and clearances.

    Coordinates are normalized by b. The initial state file is a runtime input,
    with its SHA-256 retained in the common geometry preflight. No optimization,
    change of equilibrium angles, or change of spring rest lengths occurs here.
    """
    if not cfg.flagella.initial_bundle_state_path:
        raise ValueError("phase_aligned_rear requires initial_bundle_state_path")
    path = Path(cfg.flagella.initial_bundle_state_path)
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if cfg.flagella.initial_bundle_state_sha256 != digest:
        raise ValueError(
            "Initial bundle state SHA-256 does not match pinned configuration"
        )
    data = json.loads(raw)
    if data.get("version") != 1 or data.get("units") != "b":
        raise ValueError("Unsupported initial bundle state contract")
    b_um = cfg.scale.b_um
    body = np.asarray(body_um) / b_um
    expected_body = np.asarray(data["body_positions_over_b"], dtype=float)
    if expected_body.shape != body.shape or not np.allclose(
        body, expected_body, rtol=0, atol=1e-10
    ):
        raise ValueError("Initial bundle body geometry does not match profile")
    for key, actual in (
        ("normal_bend_deg", cfg.potentials.bend.theta0_deg["normal"]),
        ("normal_torsion_deg", cfg.potentials.torsion.phi0_deg["normal"]),
        ("bead_diameter_over_b", 2 * cfg.bead_radius_m / cfg.b_m),
    ):
        if not np.isclose(data[key], actual, rtol=0, atol=1e-10):
            raise ValueError(f"Initial bundle {key} does not match profile")
    slots = list(cfg.flagella.attachment_slots or ())
    matches = [
        record for record in data["conditions"] if record["attachment_slots"] == slots
    ]
    if len(matches) != 1:
        raise ValueError(
            "Initial bundle must contain exactly one matching attachment pattern"
        )
    record = matches[0]
    points = np.asarray(record["flagella_positions_over_b"], dtype=float)
    originals = np.asarray(original_chains_um) / b_um
    if points.shape != originals.shape or not np.isfinite(points).all():
        raise ValueError("Initial bundle coordinate shape/finite check failed")
    old_bonds = np.linalg.norm(np.diff(originals, axis=1), axis=-1)
    bonds = np.linalg.norm(np.diff(points, axis=1), axis=-1)
    hook_lengths = np.linalg.norm(points[:, 0] - body[attach_ids], axis=-1)
    if not np.allclose(bonds, old_bonds, rtol=1e-8, atol=1e-10) or not np.allclose(
        hook_lengths, cfg.hook.length_over_b, rtol=1e-8, atol=1e-10
    ):
        raise ValueError(
            "Initial bundle must preserve nominal hook and shaft bond lengths"
        )
    start = int(record["distal_start_bead"])
    if not 1 <= start <= points.shape[1] - 5:
        raise ValueError("Initial bundle requires at least five distal beads")
    distal = points[:, start:]
    centered = distal - distal.mean(axis=1, keepdims=True)
    phase_error = float(np.max(np.linalg.norm(centered - centered[0], axis=-1)))
    if phase_error > 1e-9:
        raise ValueError(
            "Initial bundle distal helices do not share a common shape and phase"
        )
    tangents = np.diff(distal, axis=1)
    axes = np.cross(tangents[:, 1] - tangents[:, 0], tangents[:, 2] - tangents[:, 1])
    lengths = np.linalg.norm(axes, axis=1)
    if np.any(lengths < 1e-12):
        raise ValueError("Initial bundle distal screw axis is undefined")
    axes /= lengths[:, None]
    axes *= np.sign(np.sum(axes * (distal[:, -1] - distal[:, 0]), axis=1))[:, None]
    diameter = 2 * cfg.bead_radius_m / cfg.b_m
    min_body_bead = float(
        np.linalg.norm(points[:, :, None] - body[None, None], axis=-1).min()
    )
    min_flag_bead = np.inf
    min_flag_segment = np.inf
    for i, j in combinations(range(len(points)), 2):
        min_flag_bead = min(
            min_flag_bead,
            float(np.linalg.norm(points[i, :, None] - points[j, None], axis=-1).min()),
        )
        first = np.vstack((body[attach_ids[i]], points[i]))
        second = np.vstack((body[attach_ids[j]], points[j]))
        a, c = np.meshgrid(
            np.arange(len(first) - 1), np.arange(len(second) - 1), indexing="ij"
        )
        a, c = a.ravel(), c.ravel()
        min_flag_segment = min(
            min_flag_segment,
            float(
                np.sqrt(
                    segment_distance_squared(
                        first[a], first[a + 1], second[c], second[c + 1]
                    )
                ).min()
            ),
        )
    min_body_segment = np.inf
    for chain, attachment in zip(points, attach_ids):
        polyline = np.vstack((body[attachment], chain))
        for edge, (p, q) in enumerate(zip(polyline[:-1], polyline[1:])):
            selected = body_edges
            if edge == 0:
                selected = selected[~np.any(selected == attachment, axis=1)]
            if len(selected):
                repeated_p = np.broadcast_to(p, (len(selected), 3))
                repeated_q = np.broadcast_to(q, (len(selected), 3))
                distance = np.sqrt(
                    segment_distance_squared(
                        repeated_p,
                        repeated_q,
                        body[selected[:, 0]],
                        body[selected[:, 1]],
                    )
                )
                min_body_segment = min(min_body_segment, float(distance.min()))
    if (
        min(min_body_bead, min_flag_bead, min_body_segment, min_flag_segment)
        < diameter - 1e-8
    ):
        raise ValueError(
            "Initial bundle contains overlapping nonbonded beads or segments"
        )
    metadata = {
        "state_path": str(path),
        "state_sha256": digest,
        "distal_start_bead": start,
        "distal_bead_count": int(points.shape[1] - start),
        "common_phase_deg": record["common_phase_deg"],
        "distal_phase_shape_error_over_b": phase_error,
        "distal_screw_axis_from_rear_deg": float(
            np.degrees(np.arccos(np.clip(axes[0] @ [-1, 0, 0], -1, 1)))
        ),
        "distal_axis_unit": axes[0].tolist(),
        "min_body_bead_distance_um": min_body_bead * b_um,
        "min_flag_bead_distance_um": None
        if np.isinf(min_flag_bead)
        else min_flag_bead * b_um,
        "min_body_segment_distance_um": min_body_segment * b_um,
        "min_flag_segment_distance_um": None
        if np.isinf(min_flag_segment)
        else min_flag_segment * b_um,
        "distal_center_spread_um": float(
            np.sqrt(
                np.mean(
                    np.sum(
                        (distal.mean(axis=1) - distal.mean(axis=(0, 1))) ** 2, axis=-1
                    )
                )
            )
            * b_um
        ),
        "design_energy_over_T": record["design_energy_over_T"],
    }
    return list(points * b_um), metadata


def initial_bundle_force_metrics(cfg, model):
    """Report initial elastic stress using the production force evaluator at dt=0."""
    from sim_swim.dynamics.engine import DynamicsEngine
    from sim_swim.sim.debug_summary import _torsion_angles_rad, _triplet_angles_rad

    evaluation = DynamicsEngine(model, cfg).step(0.0)
    elastic = (
        evaluation.spring_forces
        + evaluation.bend_forces
        + evaluation.torsion_forces
        + evaluation.hook_forces
    )
    bends = model.bending_triplets[model.bending_flag_ids >= 0]
    torsions = model.torsion_quads[model.torsion_flag_ids >= 0]
    bend_error = (
        np.degrees(_triplet_angles_rad(model.positions_m, bends))
        - cfg.potentials.bend.theta0_deg["normal"]
    )
    phi = np.degrees(_torsion_angles_rad(model.positions_m, torsions))
    torsion_error = (phi - cfg.potentials.torsion.phi0_deg["normal"] + 180) % 360 - 180
    return {
        "initial_elastic_force_norm_N": float(np.linalg.norm(elastic)),
        "initial_elastic_force_max_N": float(np.linalg.norm(elastic, axis=1).max()),
        "initial_repulsion_force_norm_N": float(
            np.linalg.norm(evaluation.repulsion_forces)
        ),
        "initial_flag_bend_error_max_deg": float(abs(bend_error).max()),
        "initial_flag_torsion_error_max_deg": float(abs(torsion_error).max()),
        "initial_force_evaluation_dt_s": 0.0,
    }
