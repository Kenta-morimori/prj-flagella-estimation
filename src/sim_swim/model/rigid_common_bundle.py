"""Shape-preserving, common-screw-axis initial bundle for the 2010 hex model."""

from __future__ import annotations

import hashlib
import json
import math
from itertools import combinations
from pathlib import Path

import numpy as np

from sim_swim.model.initial_bundle import segment_distance_squared

PATTERNS = (
    (0,),
    (0, 1),
    (0, 2),
    (0, 3),
    (0, 1, 2),
    (0, 1, 3),
    (0, 1, 4),
    (0, 2, 4),
    (0, 1, 2, 3),
    (0, 1, 2, 4),
    (0, 1, 3, 4),
    (0, 1, 2, 3, 4),
    (0, 1, 2, 3, 4, 5),
)
MODE = "rigid_common_axis_rear"
PERTURBATION = "hook_axis_tilt_1deg_flag0"


def _screw_frame(points: np.ndarray) -> tuple[np.ndarray, float]:
    """Return the actual discrete helix screw axis and bead-0 radial direction."""
    bonds = np.diff(points, axis=0)
    axis = np.cross(bonds[1] - bonds[0], bonds[2] - bonds[1])
    axis /= np.linalg.norm(axis)
    axis *= np.sign(axis @ (points[-1] - points[0]))
    first = np.array([0.0, 1.0, 0.0])
    first -= axis * (first @ axis)
    first /= np.linalg.norm(first)
    second = np.cross(axis, first)
    projected = np.stack((points @ first, points @ second), axis=1)
    fit = np.linalg.lstsq(
        np.column_stack((2 * projected, np.ones(len(projected)))),
        np.sum(projected**2, axis=1),
        rcond=None,
    )[0]
    radial_components = projected[0] - fit[:2]
    radius = float(np.linalg.norm(radial_components))
    if not np.isfinite(radius) or radius <= 0:
        raise ValueError("Nominal normal helix has no finite screw radius")
    radial = (radial_components[0] * first + radial_components[1] * second) / radius
    return np.column_stack((axis, radial, np.cross(axis, radial))), radius


def _rotate_about_axis(
    values: np.ndarray, axis: np.ndarray, angle: float
) -> np.ndarray:
    return (
        values * math.cos(angle)
        + np.cross(axis, values) * math.sin(angle)
        + np.outer(values @ axis, axis) * (1 - math.cos(angle))
    )


def construct_chains(
    originals_over_b: np.ndarray,
    attachments_over_b: np.ndarray,
    slots: tuple[int, ...],
    *,
    variant: str = "base",
    hook_length_over_b: float = 0.25,
) -> tuple[np.ndarray, dict[str, float]]:
    """Move each entire nominal chain rigidly; never optimize internal beads."""
    if slots not in PATTERNS or variant not in {"base", PERTURBATION}:
        raise ValueError("Unsupported rigid common-axis attachment pattern or variant")
    if variant != "base" and slots != (0, 2, 4):
        raise ValueError("The orientation perturbation is defined only for slots024")
    if originals_over_b.shape != (len(slots), 11, 3):
        raise ValueError("Rigid common-axis bundle requires 11 beads per flagellum")
    axis_target = np.array([-1.0, 0.0, 0.0])
    result = []
    radii = []
    deltas = []
    for i, slot in enumerate(slots):
        original = originals_over_b[i]
        source_frame, radius = _screw_frame(original)
        attach = attachments_over_b[i]
        if not np.isclose(
            np.linalg.norm(attach[1:]), 0.5, atol=1e-10
        ) or not np.isclose(attach[0], 0.0, atol=1e-10):
            raise ValueError("Rigid common-axis bundle requires the central hex ring")
        ring_angle = math.atan2(attach[2], attach[1])
        cosine = (0.5**2 + radius**2 - hook_length_over_b**2) / radius
        if not -1 <= cosine <= 1:
            raise ValueError("Nominal helix radius cannot satisfy the hook length")
        delta = math.acos(cosine)
        phase = ring_angle + delta
        radial_target = np.array([0.0, math.cos(phase), math.sin(phase)])
        target_frame = np.column_stack(
            (axis_target, radial_target, np.cross(axis_target, radial_target))
        )
        rotation = target_frame @ source_frame.T
        bead0 = radius * radial_target
        chain = (original - original[0]) @ rotation.T + bead0
        if variant == PERTURBATION and i == 0:
            hook_axis = (bead0 - attach) / hook_length_over_b
            chain = bead0 + _rotate_about_axis(
                chain - bead0, hook_axis, math.radians(1.0)
            )
        result.append(chain)
        radii.append(radius)
        deltas.append(math.degrees(delta))
    return np.asarray(result), {
        "screw_radius_over_b": float(np.mean(radii)),
        "hook_azimuth_offset_deg": float(np.mean(deltas)),
        "tilt_deg": 1.0 if variant == PERTURBATION else 0.0,
    }


def _minimum_segment_gap(chains: np.ndarray, attachments: np.ndarray) -> float | None:
    minimum = math.inf
    for i, j in combinations(range(len(chains)), 2):
        left = np.vstack((attachments[i], chains[i]))
        right = np.vstack((attachments[j], chains[j]))
        for p, q in zip(left[:-1], left[1:]):
            for r, s in zip(right[:-1], right[1:]):
                minimum = min(
                    minimum,
                    math.sqrt(
                        float(
                            segment_distance_squared(
                                p[None], q[None], r[None], s[None]
                            )[0]
                        )
                    ),
                )
    return None if math.isinf(minimum) else minimum


def load_rigid_common_bundle(cfg, body_um, original_chains_um, attach_ids, body_edges):
    """Verify the pinned recipe, complete-chain shape, hook and clearances."""
    if not cfg.flagella.initial_bundle_state_path:
        raise ValueError("Rigid common-axis state path is required")
    path = Path(cfg.flagella.initial_bundle_state_path)
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != cfg.flagella.initial_bundle_state_sha256:
        raise ValueError("Rigid common-axis state SHA-256 does not match")
    data = json.loads(raw)
    if data.get("version") != 1 or data.get("kind") != MODE or data.get("units") != "b":
        raise ValueError("Unsupported rigid common-axis state contract")
    b_um = cfg.scale.b_um
    body = np.asarray(body_um) / b_um
    if not np.allclose(body, data["body_positions_over_b"], rtol=0, atol=1e-10):
        raise ValueError("Rigid common-axis body geometry does not match")
    for key, actual in (
        ("normal_bend_deg", cfg.potentials.bend.theta0_deg["normal"]),
        ("normal_torsion_deg", cfg.potentials.torsion.phi0_deg["normal"]),
        ("bead_diameter_over_b", 2 * cfg.bead_radius_m / cfg.b_m),
        ("hook_length_over_b", cfg.hook.length_over_b),
    ):
        if not np.isclose(data[key], actual, rtol=0, atol=1e-10):
            raise ValueError(f"Rigid common-axis {key} does not match profile")
    slots = tuple(cfg.flagella.attachment_slots or ())
    variant = cfg.flagella.initial_bundle_variant
    matches = [
        row
        for row in data["conditions"]
        if tuple(row["attachment_slots"]) == slots and row["variant"] == variant
    ]
    if len(matches) != 1:
        raise ValueError("Rigid common-axis state requires one matching condition")
    originals = np.asarray(original_chains_um) / b_um
    attach = body[np.asarray(attach_ids)]
    expected, recipe = construct_chains(
        originals,
        attach,
        slots,
        variant=variant,
        hook_length_over_b=cfg.hook.length_over_b,
    )
    points = np.asarray(matches[0]["flagella_positions_over_b"], dtype=float)
    if (
        points.shape != expected.shape
        or not np.isfinite(points).all()
        or not np.allclose(points, expected, rtol=0, atol=1e-9)
    ):
        raise ValueError("Rigid common-axis coordinates differ from the pinned recipe")
    pairwise_old = np.linalg.norm(originals[:, :, None] - originals[:, None], axis=-1)
    pairwise_new = np.linalg.norm(points[:, :, None] - points[:, None], axis=-1)
    if not np.allclose(pairwise_old, pairwise_new, rtol=0, atol=1e-9):
        raise ValueError("Rigid common-axis state changed the internal chain shape")
    hooks = np.linalg.norm(points[:, 0] - attach, axis=1)
    if not np.allclose(hooks, cfg.hook.length_over_b, rtol=0, atol=1e-9):
        raise ValueError("Rigid common-axis hook length changed")
    diameter = 2 * cfg.bead_radius_m / cfg.b_m
    flag_gap = _minimum_segment_gap(points, attach)
    if flag_gap is not None and flag_gap < diameter - 1e-9:
        raise ValueError("Rigid common-axis flagella segments overlap")
    body_gap = math.inf
    root_body_gap = math.inf
    self_gap = math.inf
    for chain in points:
        for i in range(len(chain) - 1):
            for j in range(i + 2, len(chain) - 1):
                self_gap = min(
                    self_gap,
                    math.sqrt(
                        float(
                            segment_distance_squared(
                                chain[i][None],
                                chain[i + 1][None],
                                chain[j][None],
                                chain[j + 1][None],
                            )[0]
                        )
                    ),
                )
    if self_gap < diameter - 1e-9:
        raise ValueError("Rigid common-axis chain self-intersects")
    # Hook and shaft segments 0..2 may traverse the represented body. Every
    # later shaft segment must clear all body spring edges by one bead diameter.
    for chain, attachment in zip(points, attach_ids):
        polyline = np.vstack((body[attachment], chain))
        for segment_index, (p, q) in enumerate(zip(polyline[:-1], polyline[1:])):
            selected = body_edges
            if segment_index == 0:
                selected = selected[~np.any(selected == attachment, axis=1)]
            distances = np.sqrt(
                segment_distance_squared(
                    np.broadcast_to(p, (len(selected), 3)),
                    np.broadcast_to(q, (len(selected), 3)),
                    body[selected[:, 0]],
                    body[selected[:, 1]],
                )
            )
            if segment_index <= 3:
                root_body_gap = min(root_body_gap, float(distances.min()))
            else:
                body_gap = min(body_gap, float(distances.min()))
    if body_gap < diameter - 1e-9:
        raise ValueError("Rigid common-axis distal chain overlaps the body")
    return list(points * b_um), {
        "state_path": str(path),
        "state_sha256": digest,
        "variant": variant,
        "whole_chain_rigid_shape_error_over_b": float(
            np.max(abs(pairwise_new - pairwise_old))
        ),
        "hook_length_error_over_b": float(np.max(abs(hooks - cfg.hook.length_over_b))),
        "min_flag_segment_distance_um": None if flag_gap is None else flag_gap * b_um,
        "min_self_segment_distance_um": self_gap * b_um,
        "min_flag_bead_distance_um": None
        if len(points) == 1
        else min(
            float(np.linalg.norm(points[i][:, None] - points[j][None], axis=-1).min())
            for i, j in combinations(range(len(points)), 2)
        )
        * b_um,
        "min_root_body_segment_distance_um": root_body_gap * b_um,
        "min_distal_body_segment_distance_um": body_gap * b_um,
        "body_penetration_allowed_through_shaft_segment": 2,
        **recipe,
    }
