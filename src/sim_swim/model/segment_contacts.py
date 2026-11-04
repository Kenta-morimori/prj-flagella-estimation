"""Shared runtime pair selection and initial segment clearance diagnostics."""

from __future__ import annotations

import numpy as np


def select_repulsion_segment_pairs(candidates, body_spring_mask, body_flagella_enabled):
    """Exclude body/body pairs; hooks are non-body segments, as in runtime."""
    first = body_spring_mask[candidates[:, 0]]
    second = body_spring_mask[candidates[:, 1]]
    bb = first & second
    bf = first ^ second
    ff = ~first & ~second
    active = ~bb
    if not body_flagella_enabled:
        active &= ~bf
    return candidates[active], {
        "body_flagella_enabled": bool(body_flagella_enabled),
        "all_non_neighbor": len(candidates),
        "body_body_excluded": int(bb.sum()),
        "body_flagella_candidate": int(bf.sum()),
        "flagella_flagella_candidate": int(ff.sum()),
        "body_flagella_active": int((active & bf).sum()),
        "flagella_flagella_active": int((active & ff).sum()),
        "active_total": int(active.sum()),
    }


def segment_clearance_metrics(positions, springs, candidates, body_mask, diameter):
    """Inspect every active nonadjacent pair using the production distance.

    Body/body and shared-endpoint exclusions come from the caller's runtime
    pair list. A witness names segment rows and bead IDs, so a root sitting on
    another body edge cannot be hidden by a bead-only clearance diagnostic.
    """
    from sim_swim.dynamics.forces import _closest_points_on_segments

    tolerance = (
        64
        * np.finfo(float).eps
        * max(diameter, float(np.ptp(positions, axis=0).max(initial=0)))
    )
    groups = {
        name: {
            "pair_count": 0,
            "overlap_count": 0,
            "min_distance_m": None,
            "witness": None,
        }
        for name in ("body_flagella", "flagella_flagella")
    }
    for i, j in candidates:
        kind = "body_flagella" if body_mask[i] ^ body_mask[j] else "flagella_flagella"
        group = groups[kind]
        distance = float(
            _closest_points_on_segments(*positions[springs[i]], *positions[springs[j]])[
                -1
            ]
        )
        if not np.isfinite(distance):
            raise ValueError("Non-finite initial segment distance")
        group["pair_count"] += 1
        group["overlap_count"] += int(distance + tolerance < diameter)
        if group["min_distance_m"] is None or distance < group["min_distance_m"]:
            group["min_distance_m"] = distance
            group["witness"] = {
                "segment_rows": [int(i), int(j)],
                "bead_indices": [springs[i].tolist(), springs[j].tolist()],
            }
    return {
        "status": "FAIL"
        if any(g["overlap_count"] for g in groups.values())
        else "PASS",
        "minimum_required_m": float(diameter),
        "roundoff_tolerance_m": tolerance,
        "groups": groups,
        "scope": "active nonadjacent spring segments, including hooks; body/body excluded",
    }


def initial_body_cylinder_clearance(body, chains, body_axis, intervals=1000):
    """Conservative all-point bound for initial filaments outside the body hull.

    A finite cylinder contains the body-centre convex hull. Its distance is
    one-Lipschitz, so subtracting half a sample spacing bounds unsampled points
    of every segment. Hooks are excluded here because they start on the body;
    their nonadjacent contacts are inspected by segment_clearance_metrics.
    This is an initial-state certificate, not a dynamic or all-spin guarantee.
    """
    axis = body_axis / np.linalg.norm(body_axis)
    center = body.mean(0)
    axial_body = (body - center) @ axis
    radial_body = body - center - np.outer(axial_body, axis)
    radius = float(np.linalg.norm(radial_body, axis=1).max())
    alpha = np.linspace(0, 1, intervals + 1)[:, None]
    bounds = []
    for points in chains:
        for a, b in zip(points[:-1], points[1:]):
            samples = (1 - alpha) * a + alpha * b - center
            axial = samples @ axis
            radial = samples - np.outer(axial, axis)
            distance = np.hypot(
                np.maximum(np.linalg.norm(radial, axis=1) - radius, 0),
                np.maximum(
                    np.maximum(axial_body.min() - axial, axial - axial_body.max()), 0
                ),
            )
            bounds.append(
                float(distance.min() - np.linalg.norm(b - a) / (2 * intervals))
            )
    return min(bounds) if bounds else None
