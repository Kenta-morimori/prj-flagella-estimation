from __future__ import annotations

import numpy as np
import pytest

from sim_swim.model.builder import ModelBuilder
from sim_swim.model.types import SimModel
from sim_swim.sim.params import SimulationConfig


def _make_cfg(
    *,
    n_prism: int = 3,
    n_flagella: int = 3,
    dz_over_b: float = 0.5,
    radius_over_b: float = 0.5,
    body_length_um: float = 2.0,
    ds_over_b: float = 0.58,
    flag_length_over_b: float = 2.32,
    seed: int = 0,
) -> SimulationConfig:
    return SimulationConfig.from_dict(
        {
            "scale": {"b_um": 1.0, "bead_radius_a_over_b": 0.1},
            "body": {
                "prism": {
                    "n_prism": n_prism,
                    "dz_over_b": dz_over_b,
                    "radius_over_b": radius_over_b,
                    "axis": "x",
                },
                "length_total_um": body_length_um,
            },
            "flagella": {
                "n_flagella": n_flagella,
                "placement_mode": "uniform",
                "discretization": {"ds_over_b": ds_over_b},
                "bond_L_over_b": ds_over_b,
                "length_over_b": flag_length_over_b,
                "helix_init": {"radius_over_b": 0.25, "pitch_over_b": 2.5},
            },
            "time": {"duration_s": 0.02, "dt_s": 1.0e-3},
            "brownian": {"enabled": False},
            "seed": {"global_seed": seed},
        }
    )


def _body_flag_pairs(model: SimModel) -> list[tuple[int, int]]:
    pairs: list[tuple[int, int]] = []
    for i, j in model.spring_pairs:
        ii = int(i)
        jj = int(j)
        i_is_body = bool(model.bead_is_body[ii])
        j_is_body = bool(model.bead_is_body[jj])
        if i_is_body ^ j_is_body:
            body_idx = ii if i_is_body else jj
            flag_idx = jj if i_is_body else ii
            pairs.append((body_idx, flag_idx))
    return pairs


def _layer_idx_for_body_bead(model: SimModel, bead_idx: int) -> int:
    for layer_idx, layer in enumerate(model.body_layer_indices):
        if int(bead_idx) in set(layer.tolist()):
            return int(layer_idx)
    raise ValueError(f"Body bead is not assigned to any layer: bead_idx={bead_idx}")


def _hook_rest_length_map(model: SimModel) -> dict[tuple[int, int], float]:
    out: dict[tuple[int, int], float] = {}
    for (i, j), rest in zip(model.spring_pairs, model.spring_rest_lengths_m):
        ii = int(i)
        jj = int(j)
        i_is_body = bool(model.bead_is_body[ii])
        j_is_body = bool(model.bead_is_body[jj])
        if not (i_is_body ^ j_is_body):
            continue
        body_idx = ii if i_is_body else jj
        flag_idx = jj if i_is_body else ii
        out[(body_idx, flag_idx)] = float(rest)
    return out


def _triplet_angle_rad(positions_m: np.ndarray, triplet: np.ndarray) -> float:
    i, j, k = (int(triplet[0]), int(triplet[1]), int(triplet[2]))
    u = positions_m[i] - positions_m[j]
    v = positions_m[k] - positions_m[j]
    nu = max(float(np.linalg.norm(u)), 1e-18)
    nv = max(float(np.linalg.norm(v)), 1e-18)
    c = float(np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0))
    return float(np.arccos(c))


def _dihedral_angle_rad(
    positions_m: np.ndarray, quad: np.ndarray | tuple[int, int, int, int]
) -> float:
    a_i, b_i, c_i, d_i = (int(quad[0]), int(quad[1]), int(quad[2]), int(quad[3]))
    a = positions_m[a_i]
    b = positions_m[b_i]
    c = positions_m[c_i]
    d = positions_m[d_i]

    b0 = a - b
    b1 = c - b
    b2 = d - c
    b1n = b1 / max(float(np.linalg.norm(b1)), 1e-18)
    v = b0 - np.dot(b0, b1n) * b1n
    w = b2 - np.dot(b2, b1n) * b1n
    x = float(np.dot(v, w))
    y = float(np.dot(np.cross(b1n, v), w))
    return float(np.arctan2(y, x))


def _wrap_deg(deg: float) -> float:
    return float((deg + 180.0) % 360.0 - 180.0)


def _estimate_helix_radius_pitch_over_b(
    model: SimModel, cfg: SimulationConfig, flag_id: int
) -> tuple[float, float]:
    flag_idx = model.flagella_indices[flag_id]
    positions_m = model.positions_m
    pts = positions_m[flag_idx]
    centered = pts - np.mean(pts, axis=0)
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    axis = vh[0]
    axis = axis / max(float(np.linalg.norm(axis)), 1e-18)
    if float(np.dot(axis, pts[-1] - pts[0])) < 0.0:
        axis = -axis

    ref = np.array([1.0, 0.0, 0.0], dtype=float)
    if abs(float(np.dot(ref, axis))) > 0.9:
        ref = np.array([0.0, 1.0, 0.0], dtype=float)
    e1 = np.cross(axis, ref)
    e1 = e1 / max(float(np.linalg.norm(e1)), 1e-18)
    e2 = np.cross(axis, e1)

    origin = pts[0]
    rel = pts - origin
    s = rel @ axis
    u = rel @ e1
    v = rel @ e2

    # Least-squares circle fit in the cross-section plane.
    mat = np.column_stack([u, v, np.ones_like(u)])
    rhs = -(u * u + v * v)
    coef, *_ = np.linalg.lstsq(mat, rhs, rcond=None)
    a, b, c = coef
    cx = -0.5 * a
    cy = -0.5 * b
    radius = np.sqrt(max(cx * cx + cy * cy - c, 0.0))

    phase = np.unwrap(np.arctan2(v - cy, u - cx))
    slope, _ = np.polyfit(phase, s, 1)
    pitch = abs(float(slope)) * 2.0 * np.pi

    return radius / cfg.b_m, pitch / cfg.b_m


def test_body_prism_bead_count_and_layout() -> None:
    cfg = _make_cfg(
        n_prism=3,
        n_flagella=2,
        dz_over_b=0.5,
        radius_over_b=0.4,
        body_length_um=2.5,
        ds_over_b=0.5,
    )

    model = ModelBuilder(cfg).build()

    assert model.body_indices.shape[0] == 3 * 6

    first_layer = model.body_layer_indices[0]
    pts = model.positions_m[first_layer] / 1e-6
    radii = np.linalg.norm(pts[:, 1:3], axis=1)
    assert np.allclose(radii, 0.4, atol=1e-6)

    x_coords = [
        np.mean(model.positions_m[layer][:, 0] / 1e-6)
        for layer in model.body_layer_indices
    ]
    diffs = np.diff(x_coords)
    assert np.allclose(diffs, 0.5, atol=1e-6)


def test_mvp_requires_n_prism_eq_3() -> None:
    cfg = _make_cfg(n_prism=4)
    with pytest.raises(ValueError, match="MVP: body.prism.n_prism must be 3"):
        ModelBuilder(cfg).build()


@pytest.mark.parametrize("n_flagella", [-1, 10])
def test_mvp_requires_n_flagella_in_range(n_flagella: int) -> None:
    cfg = _make_cfg(n_flagella=n_flagella)
    with pytest.raises(
        ValueError, match=r"MVP: flagella.n_flagella must be in \[0,9\]"
    ):
        ModelBuilder(cfg).build()


def test_body_flag_spring_connection_exists() -> None:
    cfg = _make_cfg(n_flagella=3)
    model = ModelBuilder(cfg).build()

    assert len(_body_flag_pairs(model)) == cfg.flagella.n_flagella


@pytest.mark.parametrize("n_flagella", [1, 2, 3])
def test_flagella_attach_to_center_layer_for_compat_range(n_flagella: int) -> None:
    cfg = _make_cfg(n_flagella=n_flagella)
    model = ModelBuilder(cfg).build()

    center_layer = set(
        model.body_layer_indices[len(model.body_layer_indices) // 2].tolist()
    )
    last_layer = set(model.body_layer_indices[-1].tolist())

    attached_body_indices: set[int] = set()
    for body_idx, _ in _body_flag_pairs(model):
        attached_body_indices.add(body_idx)

    assert attached_body_indices
    assert attached_body_indices.issubset(center_layer)
    assert attached_body_indices.isdisjoint(last_layer)


def test_build_supports_n_flagella_9() -> None:
    cfg = _make_cfg(n_flagella=9)
    model = ModelBuilder(cfg).build()

    assert len(model.flagella_indices) == 9


@pytest.mark.parametrize("n_flagella", [4, 5, 6, 7, 8, 9])
def test_attach_positions_have_no_duplicates_for_4_to_9(n_flagella: int) -> None:
    cfg = _make_cfg(n_flagella=n_flagella, seed=7)
    model = ModelBuilder(cfg).build()
    attached_body_indices = [body_idx for body_idx, _ in _body_flag_pairs(model)]

    assert len(attached_body_indices) == n_flagella
    assert len(set(attached_body_indices)) == n_flagella


def test_attach_positions_are_reproducible_with_same_seed() -> None:
    cfg_a = _make_cfg(n_flagella=9, seed=123)
    cfg_b = _make_cfg(n_flagella=9, seed=123)
    model_a = ModelBuilder(cfg_a).build()
    model_b = ModelBuilder(cfg_b).build()

    attach_a = sorted(body_idx for body_idx, _ in _body_flag_pairs(model_a))
    attach_b = sorted(body_idx for body_idx, _ in _body_flag_pairs(model_b))
    assert attach_a == attach_b


def test_body_flag_hook_length_is_initialized_to_0p25b() -> None:
    cfg = _make_cfg(n_flagella=3)
    model = ModelBuilder(cfg).build()

    body_flag_pairs = _body_flag_pairs(model)
    dists_m = np.asarray(
        [
            np.linalg.norm(model.positions_m[flag_idx] - model.positions_m[body_idx])
            for body_idx, flag_idx in body_flag_pairs
        ],
        dtype=float,
    )
    assert dists_m.shape[0] == cfg.flagella.n_flagella
    assert np.allclose(dists_m, 0.25 * cfg.b_m, atol=1e-12)


def test_basal_link_initial_direction_matches_local_layer_radial() -> None:
    cfg = _make_cfg(n_flagella=3)
    model = ModelBuilder(cfg).build()
    rest_map = _hook_rest_length_map(model)

    for attach_idx, first_idx, _ in model.hook_triplets:
        attach = int(attach_idx)
        first = int(first_idx)
        layer_idx = _layer_idx_for_body_bead(model, attach)
        layer = model.body_layer_indices[layer_idx]
        layer_centroid = np.mean(model.positions_m[layer], axis=0)

        link = model.positions_m[first] - model.positions_m[attach]
        link_norm = float(np.linalg.norm(link))
        assert np.isclose(link_norm, rest_map[(attach, first)], atol=1e-12)

        radial = model.positions_m[attach] - layer_centroid
        radial_unit = radial / max(float(np.linalg.norm(radial)), 1e-18)
        link_unit = link / max(link_norm, 1e-18)
        assert np.allclose(link_unit, radial_unit, atol=1e-10)


def test_spring_rest_lengths_match_initial_distances() -> None:
    cfg = _make_cfg(n_flagella=3)
    model = ModelBuilder(cfg).build()

    spring_pairs = model.spring_pairs
    actual = np.linalg.norm(
        model.positions_m[spring_pairs[:, 1]] - model.positions_m[spring_pairs[:, 0]],
        axis=1,
    )
    assert np.allclose(model.spring_rest_lengths_m, actual, atol=1e-15)


def test_flagella_intra_rest_length_is_initialized_to_0p58b() -> None:
    cfg = _make_cfg(n_flagella=2, ds_over_b=0.58)
    model = ModelBuilder(cfg).build()

    pairs = model.spring_pairs
    bi = model.bead_is_body[pairs[:, 0]]
    bj = model.bead_is_body[pairs[:, 1]]
    fi = model.bead_flag_ids[pairs[:, 0]]
    fj = model.bead_flag_ids[pairs[:, 1]]
    flag_intra_rows = np.where((~bi) & (~bj) & (fi == fj) & (fi >= 0))[0]
    assert flag_intra_rows.size > 0

    rest = model.spring_rest_lengths_m[flag_intra_rows]
    assert np.allclose(rest, 0.58 * cfg.b_m, rtol=0.0, atol=1e-12)


def test_flagellum_initial_geometry_matches_paper_values() -> None:
    cfg = _make_cfg(n_flagella=1, ds_over_b=0.58, flag_length_over_b=5.8)
    model = ModelBuilder(cfg).build()
    positions_m = model.positions_m

    pairs = model.spring_pairs
    bi = model.bead_is_body[pairs[:, 0]]
    bj = model.bead_is_body[pairs[:, 1]]
    fi = model.bead_flag_ids[pairs[:, 0]]
    fj = model.bead_flag_ids[pairs[:, 1]]
    flag_intra_rows = np.where((~bi) & (~bj) & (fi == fj) & (fi >= 0))[0]
    intra_len_over_b = (
        np.linalg.norm(
            positions_m[pairs[flag_intra_rows, 1]]
            - positions_m[pairs[flag_intra_rows, 0]],
            axis=1,
        )
        / cfg.b_m
    )
    assert abs(float(np.mean(intra_len_over_b)) - 0.58) <= 0.02 * 0.58
    assert float(np.min(intra_len_over_b)) >= 0.95 * 0.58
    assert float(np.max(intra_len_over_b)) <= 1.05 * 0.58

    radius_over_b, pitch_over_b = _estimate_helix_radius_pitch_over_b(model, cfg, 0)
    assert abs(radius_over_b - 0.25) <= 0.05 * 0.25
    assert abs(pitch_over_b - 2.5) <= 0.05 * 2.5

    bend_rows = np.where(model.bending_flag_ids >= 0)[0]
    bend_deg = np.asarray(
        [
            np.rad2deg(_triplet_angle_rad(positions_m, model.bending_triplets[row]))
            for row in bend_rows
        ],
        dtype=float,
    )
    assert abs(float(np.mean(bend_deg)) - 142.0) <= 5.0

    torsion_rows = np.where(model.torsion_flag_ids >= 0)[0]
    torsion_rad = np.asarray(
        [
            _dihedral_angle_rad(positions_m, model.torsion_quads[row])
            for row in torsion_rows
        ],
        dtype=float,
    )
    torsion_mean_deg = np.rad2deg(
        np.arctan2(np.mean(np.sin(torsion_rad)), np.mean(np.cos(torsion_rad)))
    )
    assert abs(_wrap_deg(float(torsion_mean_deg) - (-60.0))) <= 5.0


@pytest.mark.parametrize("n_flagella", [1, 3])
def test_flagella_base_tangent_points_rear(n_flagella: int) -> None:
    cfg = _make_cfg(n_flagella=n_flagella, ds_over_b=0.58)
    model = ModelBuilder(cfg).build()
    positions_m = model.positions_m

    body_points = positions_m[model.body_indices]
    centered = body_points - np.mean(body_points, axis=0)
    _, _, vh = np.linalg.svd(centered, full_matrices=False)
    body_axis = vh[0]
    body_axis = body_axis / max(float(np.linalg.norm(body_axis)), 1e-18)
    rear_dir = -body_axis

    for idx in model.flagella_indices:
        flag0 = int(idx[0])
        flag1 = int(idx[1])
        tangent0 = positions_m[flag1] - positions_m[flag0]
        tangent0 = tangent0 / max(float(np.linalg.norm(tangent0)), 1e-18)
        angle_deg = np.rad2deg(
            np.arccos(float(np.clip(np.dot(tangent0, rear_dir), -1.0, 1.0)))
        )
        assert angle_deg <= 10.0


def test_body_has_no_torsion_quads() -> None:
    cfg = _make_cfg(n_flagella=3)
    model = ModelBuilder(cfg).build()

    assert np.all(model.torsion_flag_ids >= 0)
