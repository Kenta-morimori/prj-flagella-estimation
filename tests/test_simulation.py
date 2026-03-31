from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig


def _triplet_angle_rad(positions_m: np.ndarray, triplet: np.ndarray) -> float:
    i, j, k = (int(triplet[0]), int(triplet[1]), int(triplet[2]))
    u = positions_m[i] - positions_m[j]
    v = positions_m[k] - positions_m[j]
    nu = max(float(np.linalg.norm(u)), 1e-18)
    nv = max(float(np.linalg.norm(v)), 1e-18)
    c = float(np.clip(np.dot(u, v) / (nu * nv), -1.0, 1.0))
    return float(np.arccos(c))


def _body_spring_rows(sim: Simulator) -> np.ndarray:
    pairs = sim.model.spring_pairs
    bi = sim.model.bead_is_body[pairs[:, 0]]
    bj = sim.model.bead_is_body[pairs[:, 1]]
    return np.where(bi & bj)[0]


def _body_bending_rows(sim: Simulator) -> np.ndarray:
    return np.where(sim.model.bending_flag_ids < 0)[0]


def _body_layer_lookup(sim: Simulator) -> dict[int, np.ndarray]:
    out: dict[int, np.ndarray] = {}
    for layer in sim.model.body_layer_indices:
        layer_int = layer.astype(int, copy=False)
        for bead_idx in layer_int:
            out[int(bead_idx)] = layer_int
    return out


def _hook_rest_lookup(sim: Simulator) -> dict[tuple[int, int], float]:
    out: dict[tuple[int, int], float] = {}
    for (i, j), rest in zip(sim.model.spring_pairs, sim.model.spring_rest_lengths_m):
        ii = int(i)
        jj = int(j)
        i_is_body = bool(sim.model.bead_is_body[ii])
        j_is_body = bool(sim.model.bead_is_body[jj])
        if not (i_is_body ^ j_is_body):
            continue
        body_idx = ii if i_is_body else jj
        flag_idx = jj if i_is_body else ii
        out[(body_idx, flag_idx)] = float(rest)
    return out


def _assert_basal_link_geometry(sim: Simulator, positions_m: np.ndarray) -> None:
    layer_lookup = _body_layer_lookup(sim)
    rest_lookup = _hook_rest_lookup(sim)
    for attach_idx, first_idx, _ in sim.model.hook_triplets:
        attach = int(attach_idx)
        first = int(first_idx)
        layer = layer_lookup[attach]
        layer_centroid = np.mean(positions_m[layer], axis=0)

        link = positions_m[first] - positions_m[attach]
        link_norm = float(np.linalg.norm(link))
        expected_rest = rest_lookup[(attach, first)]
        assert np.isclose(link_norm, expected_rest, atol=1e-12)

        radial = positions_m[attach] - layer_centroid
        radial_unit = radial / max(float(np.linalg.norm(radial)), 1e-18)
        link_unit = link / max(link_norm, 1e-18)
        assert np.allclose(link_unit, radial_unit, atol=1e-10)


def _measure_body_constraint_drift(sim: Simulator, n_steps: int) -> tuple[float, float]:
    model = sim.model
    body_spring_rows = _body_spring_rows(sim)
    body_bending_rows = _body_bending_rows(sim)
    theta0, _ = sim.engine._state_angles_rad()
    theta0_body = theta0[body_bending_rows].copy()
    spring_l0 = model.spring_rest_lengths_m[body_spring_rows].copy()
    spring_pairs = model.spring_pairs[body_spring_rows]

    max_stretch_ratio = 0.0
    max_angle_error_rad = 0.0
    for _ in range(max(1, int(n_steps))):
        sim.engine.step(sim.config.dt_star)
        pos = model.positions_m

        spring_l = np.linalg.norm(
            pos[spring_pairs[:, 1]] - pos[spring_pairs[:, 0]],
            axis=1,
        )
        stretch_ratio = np.abs(spring_l - spring_l0) / np.maximum(spring_l0, 1e-30)
        max_stretch_ratio = max(max_stretch_ratio, float(np.max(stretch_ratio)))

        theta = np.asarray(
            [
                _triplet_angle_rad(pos, model.bending_triplets[row])
                for row in body_bending_rows
            ],
            dtype=float,
        )
        angle_error_rad = np.abs(theta - theta0_body)
        max_angle_error_rad = max(max_angle_error_rad, float(np.max(angle_error_rad)))

    return max_stretch_ratio, max_angle_error_rad


def _make_cfg(
    motor_torque_Nm: float = -1.0,
    hook_enabled: bool = True,
    n_flagella: int = 3,
) -> SimulationConfig:
    return SimulationConfig.from_dict(
        {
            "scale": {"b_um": 1.0, "bead_radius_a_over_b": 0.1},
            "body": {
                "prism": {
                    "n_prism": 3,
                    "dz_over_b": 0.5,
                    "radius_over_b": 0.5,
                    "axis": "x",
                },
                "length_total_um": 2.0,
            },
            "flagella": {
                "n_flagella": n_flagella,
                "placement_mode": "uniform",
                "discretization": {"ds_over_b": 0.58},
                "bond_L_over_b": 0.58,
                "length_over_b": 2.32,
                "helix_init": {"radius_over_b": 0.25, "pitch_over_b": 2.5},
            },
            "fluid": {"viscosity_Pa_s": 1.0e-3},
            "motor": {"torque_Nm": motor_torque_Nm, "reverse_n_flagella": 1},
            "potentials": {
                "spring": {"H_over_T_over_b": 10.0, "s": 0.1},
                "bend": {
                    "kb_over_T": 20.0,
                    "theta0_deg": {
                        "normal": 142.0,
                        "semicoiled": 90.0,
                        "curly1": 105.0,
                    },
                },
                "torsion": {
                    "kt_over_T": 10.0,
                    "phi0_deg": {"normal": -60.0, "semicoiled": 65.0, "curly1": 120.0},
                },
                "spring_spring_repulsion": {
                    "A_ss_over_T": 1.0,
                    "a_ss_over_b": 0.2,
                    "cutoff_over_b": 0.2,
                },
            },
            "hook": {"enabled": hook_enabled, "threshold_deg": 90.0, "kb_over_T": 20.0},
            "run_tumble": {
                "run_tau": 20.0,
                "tumble_tau": 8.0,
                "semicoiled_tau": 4.0,
                "curly1_tau": 4.0,
            },
            "time": {"duration_s": 5.0e-5, "dt_s": 1.0e-3},
            "output_sampling": {"out_all_steps_3d": True, "fps_out_2d": 25.0},
            "brownian": {
                "enabled": False,
                "temperature_K": 298.0,
                "method": "cholesky",
                "jitter": 1.0e-20,
            },
            "render": {
                "image_size_px": 128,
                "pixel_size_um": 0.203,
                "flagella_linewidth_px": 2.0,
                "render_flagella": True,
                "save_frames_3d": False,
                "follow_camera_3d": True,
                "view_range_um": 8.0,
                "timestamp_3d": False,
                "label_flagella": True,
                "follow_camera_2d": False,
                "save_frames_2d": False,
            },
            "seed": {"global_seed": 0},
            "output": {"base_dir": "outputs"},
        }
    )


def _run_and_load_step_summary(
    sim: Simulator, duration_s: float, summary_dir: Path
) -> list[dict[str, str]]:
    sim.run(duration_s, step_summary_dir=summary_dir)
    csv_path = summary_dir / "step_summary.csv"
    assert csv_path.is_file()
    with csv_path.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    assert rows
    return rows


def _assert_hook_stats_over_all_steps(
    rows: list[dict[str, str]], b_um: float, n_flagella: int
) -> None:
    target_um = 0.25 * b_um
    mean_tol_um = 0.05 * target_um
    min_limit_um = 0.85 * target_um
    max_limit_um = 1.15 * target_um

    for row in rows:
        assert int(row["hook_count"]) == n_flagella
        hook_mean = float(row["hook_len_mean_um"])
        hook_min = float(row["hook_len_min_um"])
        hook_max = float(row["hook_len_max_um"])
        assert abs(hook_mean - target_um) <= mean_tol_um
        assert hook_min >= min_limit_um
        assert hook_max <= max_limit_um


def _assert_flag_intra_stats_over_all_steps(
    rows: list[dict[str, str]], b_um: float
) -> None:
    target_um = 0.58 * b_um
    mean_tol_um = 0.02 * target_um
    min_limit_um = 0.95 * target_um
    max_limit_um = 1.05 * target_um

    for row in rows:
        assert int(row["flag_intra_count"]) > 0
        intra_mean = float(row["flag_intra_len_mean_um"])
        intra_min = float(row["flag_intra_len_min_um"])
        intra_max = float(row["flag_intra_len_max_um"])
        assert abs(intra_mean - target_um) <= mean_tol_um
        assert intra_min >= min_limit_um
        assert intra_max <= max_limit_um


def _projection_all_off_overrides() -> dict[str, dict[str, bool]]:
    return {
        "projection": {
            "enable_body_rigid_projection": False,
            "enable_hook_length_projection": False,
            "enable_basal_link_direction_projection": False,
            "enable_flagella_chain_length_projection": False,
            "enable_flagella_template_projection": False,
        }
    }


def _assert_flag_angle_error_over_all_steps(
    rows: list[dict[str, str]],
    bend_mean_deg_limit: float,
    bend_max_deg_limit: float,
    torsion_mean_deg_limit: float,
    torsion_max_deg_limit: float,
) -> None:
    for row in rows:
        bend_mean = float(row["flag_bend_err_mean_deg"])
        bend_max = float(row["flag_bend_err_max_deg"])
        torsion_mean = float(row["flag_torsion_err_mean_deg"])
        torsion_max = float(row["flag_torsion_err_max_deg"])
        assert bend_mean <= bend_mean_deg_limit
        assert bend_max <= bend_max_deg_limit
        assert torsion_mean <= torsion_mean_deg_limit
        assert torsion_max <= torsion_max_deg_limit


def test_short_run_no_nan_inf() -> None:
    cfg = _make_cfg()
    sim = Simulator(cfg)
    states = sim.run(cfg.time.duration_s)

    assert len(states) >= 2

    arr_pos = np.array([s.position_um for s in states], dtype=float)
    arr_q = np.array([s.quaternion for s in states], dtype=float)
    arr_v = np.array([s.velocity_um_s for s in states], dtype=float)
    arr_w = np.array([s.omega_rad_s for s in states], dtype=float)

    assert np.isfinite(arr_pos).all()
    assert np.isfinite(arr_q).all()
    assert np.isfinite(arr_v).all()
    assert np.isfinite(arr_w).all()


def test_run_writes_step_summary_csv(tmp_path: Path) -> None:
    cfg = _make_cfg()
    sim = Simulator(cfg)
    sim.run(cfg.time.duration_s, step_summary_dir=tmp_path / "sim")

    step_csv = tmp_path / "sim" / "step_summary.csv"
    step_full_csv = tmp_path / "sim" / "step_summary_full.csv"
    assert step_csv.is_file()
    assert not step_full_csv.exists()

    with step_csv.open("r", encoding="utf-8", newline="") as f:
        rows = list(csv.DictReader(f))
    assert rows
    first = rows[0]
    assert "step" in first
    assert "F_total_mean_all" in first
    assert "hook_count" in first
    assert "flag_intra_count" in first
    assert "flag_bend_err_mean_deg" in first
    assert "flag_torsion_err_mean_deg" in first
    assert "flag_state_changed" in first
    assert "projection_body_rigid_enabled" in first
    assert "projection_hook_length_enabled" in first
    assert "projection_basal_link_direction_enabled" in first
    assert "projection_flagella_chain_length_enabled" in first
    assert "projection_flagella_template_enabled" in first
    assert "torsion_fd_eps_m" in first
    assert "torsion_fd_eps_over_b" in first
    assert "hook_len_mean_over_b" in first
    assert "hook_len_rel_err_mean" in first
    assert "hook_angle_mean_deg" in first
    assert "hook_angle_err_max_deg" in first
    assert "flag_bond_len_mean_over_b" in first
    assert "flag_bond_rel_err_mean" in first
    assert "motor_ra_len_um" in first
    assert "motor_rb_len_um" in first
    assert "motor_Ta_norm" in first
    assert "motor_Tb_norm" in first
    assert "motor_Fa_norm" in first
    assert "motor_Fb_norm" in first
    assert "motor_axis_vs_rear_direction_angle_deg" in first
    assert "motor_attach_force_norm" in first
    assert "motor_first_force_norm" in first
    assert "motor_second_force_norm" in first
    assert first["brownian_enabled"] in {"False", "false", "0"}


def test_run_writes_initial_geometry_summary_json(tmp_path: Path) -> None:
    cfg = _make_cfg(
        motor_torque_Nm=0.0, hook_enabled=True, n_flagella=1
    ).with_overrides(
        {
            "flagella": {
                "init_mode": "paper_table1",
                "n_beads_per_flagellum": 11,
            }
        }
    )
    sim = Simulator(cfg)
    sim.run(cfg.time.duration_s, step_summary_dir=tmp_path / "sim")

    summary_path = tmp_path / "sim" / "initial_geometry_summary.json"
    assert summary_path.is_file()
    data = json.loads(summary_path.read_text(encoding="utf-8"))
    assert data["flagella"]["init_mode"] == "paper_table1"
    assert data["flagella"]["n_flagella"] == 1
    assert data["flagella"]["n_beads_per_flagellum"] == 11
    assert len(data["per_flagellum"]) == 1

    flag0 = data["per_flagellum"][0]
    assert flag0["bead_count"] == 11
    assert "contour_length_um" in flag0
    assert "initial_bend_err_max_deg" in flag0
    assert "initial_torsion_err_max_deg" in flag0


@pytest.mark.parametrize(
    ("disabled_key", "expected_counter_key"),
    [
        ("enable_hook_length_projection", "hook_length"),
        ("enable_basal_link_direction_projection", "basal_direction"),
        ("enable_flagella_chain_length_projection", "flag_chain"),
        ("enable_flagella_template_projection", "flag_template"),
    ],
)
def test_projection_toggle_bypasses_target_projection(
    disabled_key: str,
    expected_counter_key: str,
) -> None:
    cfg = _make_cfg(
        motor_torque_Nm=0.0, hook_enabled=True, n_flagella=1
    ).with_overrides(
        {
            "projection": {
                disabled_key: False,
            }
        }
    )
    sim = Simulator(cfg)

    counters = {
        "hook_length": 0,
        "basal_direction": 0,
        "flag_chain": 0,
        "flag_template": 0,
    }

    original_distance_pairs = sim.engine._project_distance_pairs
    original_basal = sim.engine._project_basal_link_direction
    original_chain = sim.engine._project_flagella_chain_lengths
    original_template = sim.engine._project_flagella_template

    def wrapped_distance_pairs(
        positions_m: np.ndarray,
        pair_rows: np.ndarray,
        iterations: int,
        fixed_mask: np.ndarray | None = None,
    ) -> np.ndarray:
        if np.array_equal(pair_rows, sim.engine.hook_spring_rows):
            counters["hook_length"] += 1
        return original_distance_pairs(positions_m, pair_rows, iterations, fixed_mask)

    def wrapped_basal(positions_m: np.ndarray) -> np.ndarray:
        counters["basal_direction"] += 1
        return original_basal(positions_m)

    def wrapped_chain(positions_m: np.ndarray, iterations: int) -> np.ndarray:
        counters["flag_chain"] += 1
        return original_chain(positions_m, iterations)

    def wrapped_template(positions_m: np.ndarray) -> np.ndarray:
        counters["flag_template"] += 1
        return original_template(positions_m)

    sim.engine._project_distance_pairs = wrapped_distance_pairs
    sim.engine._project_basal_link_direction = wrapped_basal
    sim.engine._project_flagella_chain_lengths = wrapped_chain
    sim.engine._project_flagella_template = wrapped_template

    sim.engine.step(cfg.dt_star)

    assert counters[expected_counter_key] == 0
    for key, val in counters.items():
        if key == expected_counter_key:
            continue
        assert val > 0


@pytest.mark.parametrize(
    ("motor_torque", "hook_enabled"),
    [
        (0.0, False),
        (0.0, True),
        (-1.0, True),
    ],
)
def test_body_flag_bond_metric_is_numeric_across_modes(
    tmp_path: Path, motor_torque: float, hook_enabled: bool
) -> None:
    cfg = _make_cfg(motor_torque_Nm=motor_torque, hook_enabled=hook_enabled)
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, cfg.time.duration_s, tmp_path / "sim")

    first = rows[0]
    assert first["pos_all_finite"] in {"True", "true", "1"}
    hook_len = float(first["hook_len_mean_um"])
    assert np.isfinite(hook_len)


@pytest.mark.parametrize("hook_enabled", [False, True])
def test_zero_torque_initial_steps_do_not_scatter(hook_enabled: bool) -> None:
    cfg = _make_cfg(motor_torque_Nm=0.0, hook_enabled=hook_enabled)
    sim = Simulator(cfg)

    mean_disp_um: list[float] = []
    for _ in range(5):
        diag = sim.engine.step(cfg.dt_star)
        disp = np.linalg.norm(diag.positions_after_m - diag.positions_before_m, axis=1)
        mean_disp_um.append(float(np.mean(disp) * 1e6))

    assert max(mean_disp_um) < 0.7


def test_hook_bond_stats_columns_and_expected_values(tmp_path: Path) -> None:
    cfg = _make_cfg(motor_torque_Nm=0.0, hook_enabled=True, n_flagella=3)
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, cfg.time.duration_s, tmp_path / "sim")

    first = rows[0]
    assert int(first["hook_count"]) == cfg.flagella.n_flagella
    mean_len = float(first["hook_len_mean_um"])
    min_len = float(first["hook_len_min_um"])
    max_len = float(first["hook_len_max_um"])
    target_um = 0.25 * cfg.scale.b_um
    tol_um = 0.05 * target_um
    assert abs(mean_len - target_um) <= tol_um
    assert abs(min_len - target_um) <= tol_um
    assert abs(max_len - target_um) <= tol_um


def test_basal_link_initial_geometry_matches_local_radial() -> None:
    cfg = _make_cfg(motor_torque_Nm=0.0, hook_enabled=True, n_flagella=3)
    sim = Simulator(cfg)
    _assert_basal_link_geometry(sim, sim.model.positions_m)


def test_basal_link_geometry_after_one_step() -> None:
    cfg = _make_cfg(motor_torque_Nm=1.0e-18, hook_enabled=True, n_flagella=3)
    sim = Simulator(cfg)
    sim.engine.step(cfg.dt_star)
    _assert_basal_link_geometry(sim, sim.model.positions_m)


def test_basal_link_geometry_after_multi_step() -> None:
    cfg = _make_cfg(motor_torque_Nm=1.0e-18, hook_enabled=True, n_flagella=3)
    sim = Simulator(cfg)

    for _ in range(300):
        sim.engine.step(cfg.dt_star)

    _assert_basal_link_geometry(sim, sim.model.positions_m)


def test_basal_projection_directly_moves_first_bead_only() -> None:
    cfg = _make_cfg(motor_torque_Nm=1.0e-18, hook_enabled=True, n_flagella=3)
    sim = Simulator(cfg)
    before = sim.model.positions_m.copy()
    perturbed = before.copy()

    all_flag_idx = np.concatenate(sim.model.flagella_indices).astype(int, copy=False)
    for bead_idx in all_flag_idx:
        perturbed[bead_idx] += np.array([1.0e-8, -2.0e-8, 1.5e-8], dtype=float)

    projected = sim.engine._project_basal_link_direction(perturbed)

    first_idx = {int(t[1]) for t in sim.model.hook_triplets.tolist()}
    non_first_flag_idx = sorted(set(all_flag_idx.tolist()) - first_idx)
    if non_first_flag_idx:
        assert np.allclose(
            projected[non_first_flag_idx],
            perturbed[non_first_flag_idx],
            atol=0.0,
            rtol=0.0,
        )


def test_body_constraints_match_initial_reference() -> None:
    cfg = _make_cfg(motor_torque_Nm=0.0, hook_enabled=False, n_flagella=0)
    sim = Simulator(cfg)
    model = sim.model

    body_spring_rows = _body_spring_rows(sim)
    spring_pairs = model.spring_pairs[body_spring_rows]
    spring_l = np.linalg.norm(
        model.positions_m[spring_pairs[:, 1]] - model.positions_m[spring_pairs[:, 0]],
        axis=1,
    )
    spring_l0 = model.spring_rest_lengths_m[body_spring_rows]
    rel_err = np.max(np.abs(spring_l - spring_l0) / np.maximum(spring_l0, 1e-30))
    assert rel_err < 1e-12

    body_bending_rows = _body_bending_rows(sim)
    theta = np.asarray(
        [
            _triplet_angle_rad(model.positions_m, model.bending_triplets[row])
            for row in body_bending_rows
        ],
        dtype=float,
    )
    theta0, _ = sim.engine._state_angles_rad()
    theta0_body = theta0[body_bending_rows]
    max_theta_err = float(np.max(np.abs(theta - theta0_body)))
    assert max_theta_err < 1e-12


def test_body_constraints_nearly_invariant_after_one_step_without_external_forces() -> (
    None
):
    cfg = _make_cfg(motor_torque_Nm=0.0, hook_enabled=False, n_flagella=0)
    sim = Simulator(cfg)
    model = sim.model

    body_spring_rows = _body_spring_rows(sim)
    body_bending_rows = _body_bending_rows(sim)
    theta0, _ = sim.engine._state_angles_rad()
    theta0_body = theta0[body_bending_rows].copy()
    spring_l0 = model.spring_rest_lengths_m[body_spring_rows].copy()

    sim.engine.step(cfg.dt_star)
    pos = model.positions_m

    spring_pairs = model.spring_pairs[body_spring_rows]
    spring_l = np.linalg.norm(pos[spring_pairs[:, 1]] - pos[spring_pairs[:, 0]], axis=1)
    rel_err = np.max(np.abs(spring_l - spring_l0) / np.maximum(spring_l0, 1e-30))
    assert rel_err < 1e-4

    theta = np.asarray(
        [
            _triplet_angle_rad(pos, model.bending_triplets[row])
            for row in body_bending_rows
        ],
        dtype=float,
    )
    max_theta_err = float(np.max(np.abs(theta - theta0_body)))
    assert max_theta_err < 1e-3


def test_body_rigidity_is_preserved_over_multiple_steps_without_external_forces() -> (
    None
):
    cfg = _make_cfg(motor_torque_Nm=0.0, hook_enabled=False, n_flagella=0)
    sim = Simulator(cfg)
    max_stretch_ratio, max_angle_error_rad = _measure_body_constraint_drift(
        sim, n_steps=1000
    )

    assert max_stretch_ratio < 1e-4, (
        f"body distance constraint drifted: max_stretch_ratio={max_stretch_ratio:.6e}"
    )
    assert max_angle_error_rad < 1e-3, (
        f"body angle constraint drifted: max_angle_error_rad={max_angle_error_rad:.6e}"
    )


def test_motor_on_1e19_stability_and_body_constraints() -> None:
    cfg = _make_cfg(motor_torque_Nm=1.0e-19, hook_enabled=True, n_flagella=1)
    sim = Simulator(cfg)
    max_stretch_ratio, max_angle_error_rad = _measure_body_constraint_drift(
        sim, n_steps=2000
    )

    assert np.isfinite(sim.model.positions_m).all()
    assert max_stretch_ratio < 5e-4, (
        f"motor-on body distance drifted: max_stretch_ratio={max_stretch_ratio:.6e}"
    )
    assert max_angle_error_rad < 2e-3, (
        f"motor-on body angle drifted: max_angle_error_rad={max_angle_error_rad:.6e}"
    )


def test_hook_length_multi_step_motor_off(tmp_path: Path) -> None:
    cfg = _make_cfg(
        motor_torque_Nm=0.0, hook_enabled=True, n_flagella=1
    ).with_overrides({"time": {"duration_s": 2.0}})
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, 2.0, tmp_path / "sim")

    assert len(rows) >= 2000
    _assert_hook_stats_over_all_steps(rows, cfg.scale.b_um, cfg.flagella.n_flagella)


def test_projection_all_off_motor_off_has_stable_short_run_metrics(
    tmp_path: Path,
) -> None:
    cfg = _make_cfg(
        motor_torque_Nm=0.0, hook_enabled=True, n_flagella=1
    ).with_overrides(
        {
            "time": {"duration_s": 0.002},
            **_projection_all_off_overrides(),
        }
    )
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, cfg.time.duration_s, tmp_path / "sim")

    assert len(rows) >= 2
    for row in rows:
        assert row["pos_all_finite"] in {"True", "true", "1"}
        assert int(row["hook_count"]) == cfg.flagella.n_flagella
        assert int(row["flag_intra_count"]) > 0

        hook_mean_over_b = float(row["hook_len_mean_over_b"])
        hook_max_over_b = float(row["hook_len_max_over_b"])
        flag_min_over_b = float(row["flag_bond_len_min_over_b"])
        flag_max_over_b = float(row["flag_bond_len_max_over_b"])

        # projectionを全OFFにした短時間での「極端崩壊なし」を確認する境界。
        assert np.isfinite(hook_mean_over_b)
        assert np.isfinite(hook_max_over_b)
        assert np.isfinite(flag_min_over_b)
        assert np.isfinite(flag_max_over_b)
        assert hook_max_over_b <= 10.0
        assert flag_max_over_b <= 10.0


def test_projection_all_off_motor_on_outputs_comparable_diagnostics(
    tmp_path: Path,
) -> None:
    cfg = _make_cfg(
        motor_torque_Nm=1.0e-18, hook_enabled=True, n_flagella=1
    ).with_overrides(
        {
            "time": {"duration_s": 0.01},
            **_projection_all_off_overrides(),
        }
    )
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, cfg.time.duration_s, tmp_path / "sim")

    assert rows
    first = rows[0]
    assert first["projection_body_rigid_enabled"] in {"False", "false", "0"}
    assert first["projection_hook_length_enabled"] in {"False", "false", "0"}
    assert first["projection_basal_link_direction_enabled"] in {
        "False",
        "false",
        "0",
    }
    assert first["projection_flagella_chain_length_enabled"] in {
        "False",
        "false",
        "0",
    }
    assert first["projection_flagella_template_enabled"] in {
        "False",
        "false",
        "0",
    }
    assert np.isfinite(float(first["hook_len_mean_over_b"]))
    assert np.isfinite(float(first["flag_bond_len_mean_over_b"]))
    assert np.isfinite(float(first["motor_ra_len_um"]))
    assert np.isfinite(float(first["motor_rb_len_um"]))
    assert np.isfinite(float(first["motor_Ta_norm"]))
    assert np.isfinite(float(first["motor_Tb_norm"]))
    assert np.isfinite(float(first["motor_Fa_norm"]))
    assert np.isfinite(float(first["motor_Fb_norm"]))
    assert np.isfinite(float(first["motor_axis_vs_rear_direction_angle_deg"]))
    assert np.isfinite(float(first["motor_attach_force_norm"]))
    assert np.isfinite(float(first["motor_first_force_norm"]))
    assert np.isfinite(float(first["motor_second_force_norm"]))


def test_torsion_fd_eps_sweep_is_traceable_and_reduces_step0_torsion_force(
    tmp_path: Path,
) -> None:
    eps_values = [0.1, 0.001, 0.0001]
    torsion_force_step0: list[float] = []

    for eps in eps_values:
        cfg = _make_cfg(
            motor_torque_Nm=0.0,
            hook_enabled=True,
            n_flagella=1,
        ).with_overrides(
            {
                "time": {"duration_s": 0.002},
                "flagella": {
                    "init_mode": "paper_table1",
                    "n_beads_per_flagellum": 11,
                },
                "potentials": {"torsion": {"fd_eps_over_b": eps}},
                **_projection_all_off_overrides(),
            }
        )
        sim = Simulator(cfg)
        rows = _run_and_load_step_summary(
            sim,
            cfg.time.duration_s,
            tmp_path / f"sim_eps_{eps}",
        )
        first = rows[0]
        assert float(first["torsion_fd_eps_over_b"]) == pytest.approx(eps)
        torsion_force_step0.append(float(first["F_torsion_mean_flag"]))

    assert torsion_force_step0[1] < torsion_force_step0[0]
    assert torsion_force_step0[2] < torsion_force_step0[0]


def test_hook_length_multi_step_motor_on(tmp_path: Path) -> None:
    cfg = _make_cfg(
        motor_torque_Nm=1.0e-18,
        hook_enabled=True,
        n_flagella=1,
    ).with_overrides({"time": {"duration_s": 2.0}})
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, 2.0, tmp_path / "sim")

    assert len(rows) >= 2000
    _assert_hook_stats_over_all_steps(rows, cfg.scale.b_um, cfg.flagella.n_flagella)


def test_flag_intra_length_multi_step_motor_off(tmp_path: Path) -> None:
    cfg = _make_cfg(
        motor_torque_Nm=0.0, hook_enabled=True, n_flagella=1
    ).with_overrides({"time": {"duration_s": 2.0}})
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, 2.0, tmp_path / "sim")

    assert len(rows) >= 2000
    _assert_flag_intra_stats_over_all_steps(rows, cfg.scale.b_um)


def test_flag_intra_length_multi_step_motor_on(tmp_path: Path) -> None:
    cfg = _make_cfg(
        motor_torque_Nm=1.0e-18,
        hook_enabled=True,
        n_flagella=1,
    ).with_overrides({"time": {"duration_s": 2.0}})
    sim = Simulator(cfg)
    rows = _run_and_load_step_summary(sim, 2.0, tmp_path / "sim")

    assert len(rows) >= 2000
    _assert_flag_intra_stats_over_all_steps(rows, cfg.scale.b_um)
    _assert_flag_angle_error_over_all_steps(
        rows,
        bend_mean_deg_limit=6.0,
        bend_max_deg_limit=20.0,
        torsion_mean_deg_limit=6.0,
        torsion_max_deg_limit=25.0,
    )
