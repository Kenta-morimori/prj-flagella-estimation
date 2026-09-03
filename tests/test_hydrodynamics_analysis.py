from pathlib import Path
import subprocess
import sys

import numpy as np

from sim_swim.analysis.hydrodynamics import (
    HYDRO_ARCHIVE_FORMAT,
    HydroSample,
    load_hydro_archive,
    rpy_flow_velocity,
    save_hydro_archive,
    stokes_fluid_resistance,
    velocity_contributions,
)
from sim_swim.dynamics.hydro_rpy import compute_rpy_pair_mobility
from sim_swim.analysis.hydrodynamics_campaign import (
    FLOW_SLICE_GRID_SIZE,
    analyze_campaign,
    body_fixed_flow_slice,
)
from sim_swim.analysis.hydrodynamics_replay import _visible_with_common_reference
from sim_swim.analysis.flagella_count_behavior import save_state_archive
from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    load_yaml,
    normalize_campaign_config,
)
from sim_swim.sim.core import Simulator
from sim_swim.sim.core import SimulationState
from sim_swim.sim.params import SimulationConfig


def test_velocity_contributions_are_additive_and_separate_other_beads() -> None:
    positions = np.asarray([[0.0, 0.0, 0.0], [3.0e-6, 0.0, 0.0]])
    forces = np.asarray([[1.0e-12, 0.0, 0.0], [0.0, 2.0e-12, 0.0]])
    result = velocity_contributions(
        positions, forces, bead_radius_m=1.0e-6, viscosity_Pa_s=1.0e-3
    )
    np.testing.assert_allclose(result["source"], result["self"] + result["other"])
    assert np.linalg.norm(result["other"]) > 0.0


def test_rpy_flow_is_linear_in_source_selection() -> None:
    positions = np.asarray([[0.0, 0.0, 0.0], [3.0e-6, 0.0, 0.0]])
    forces = np.asarray([[1.0e-12, 0.0, 0.0], [0.0, 2.0e-12, 0.0]])
    points = np.asarray([[0.0, 4.0e-6, 0.0]])
    kwargs = dict(bead_radius_m=1.0e-6, viscosity_Pa_s=1.0e-3)
    total = rpy_flow_velocity(points, positions, forces, **kwargs)
    split = sum(
        (
            rpy_flow_velocity(points, positions, forces, source_mask=mask, **kwargs)
            for mask in (np.asarray([True, False]), np.asarray([False, True]))
        ),
        start=np.zeros_like(total),
    )
    np.testing.assert_allclose(total, split)


def test_dense_flow_reconstruction_matches_pair_kernel() -> None:
    positions = np.asarray([[0.0, 0.0, 0.0], [3.0e-6, 0.0, 0.0]])
    forces = np.asarray([[1.0e-12, 0.0, 0.0], [0.0, 2.0e-12, 0.0]])
    points = np.asarray([[0.0, 4.0e-6, 0.0], [2.0e-6, -1.0e-6, 0.0]])
    kwargs = dict(bead_radius_m=1.0e-6, viscosity_Pa_s=1.0e-3)
    expected = np.asarray(
        [
            sum(
                (
                    compute_rpy_pair_mobility(point - position, **kwargs) @ force
                    for position, force in zip(positions, forces, strict=True)
                ),
                start=np.zeros(3),
            )
            for point in points
        ]
    )
    np.testing.assert_allclose(
        rpy_flow_velocity(points, positions, forces, **kwargs), expected
    )


def test_stokes_fluid_resistance_is_exact_force_balance() -> None:
    mechanical = np.asarray([[1.0e-12, -2.0e-12, 0.0], [-3.0e-12, 4.0e-12, 0.0]])
    np.testing.assert_allclose(mechanical + stokes_fluid_resistance(mechanical), 0.0)


def test_individual_source_contributions_sum_to_total_velocity() -> None:
    positions = np.asarray([[0.0, 0.0, 0.0], [3.0e-6, 0.0, 0.0], [0.0, 3.0e-6, 0.0]])
    forces = np.asarray([[1.0e-12, 0.0, 0.0], [0.0, 2.0e-12, 0.0], [0.0, 0.0, 3.0e-12]])
    kwargs = dict(bead_radius_m=1.0e-6, viscosity_Pa_s=1.0e-3)
    total = velocity_contributions(positions, forces, **kwargs)["source"]
    split = sum(
        (
            velocity_contributions(positions, forces, source_mask=mask, **kwargs)[
                "source"
            ]
            for mask in (
                np.asarray([True, False, False]),
                np.asarray([False, True, False]),
                np.asarray([False, False, True]),
            )
        ),
        start=np.zeros_like(total),
    )
    np.testing.assert_allclose(total, split)


def test_source_contributions_share_one_visual_velocity_reference() -> None:
    small = np.asarray([[1.0, 0.0, 0.0]])
    large = np.asarray([[10.0, 0.0, 0.0]])
    displayed_small, displayed_large = _visible_with_common_reference([small, large])
    np.testing.assert_allclose(
        np.linalg.norm(displayed_large), np.linalg.norm(displayed_small) * 10.0
    )


def test_hydro_archive_round_trip(tmp_path: Path) -> None:
    path = tmp_path / "hydro_archive.npz"
    save_hydro_archive(
        path,
        [
            HydroSample(0.0, np.zeros((2, 3)), np.ones((2, 3))),
            HydroSample(0.001, np.ones((2, 3)), np.zeros((2, 3))),
        ],
        bead_is_body=np.asarray([True, False]),
        bead_flagella_id=np.asarray([-1, 0]),
        bead_radius_m=1.0e-6,
        viscosity_Pa_s=1.0e-3,
        provenance={"commit": "abc"},
    )
    archive = load_hydro_archive(path)
    assert HYDRO_ARCHIVE_FORMAT == "sim_swim.hydro_archive"
    np.testing.assert_allclose(archive.t_s, [0.0, 0.001])
    assert archive.provenance == {"commit": "abc"}


def test_opt_in_simulator_records_compact_force_position_samples() -> None:
    raw = load_yaml(Path("conf/sim_swim_2010.yaml"))
    cfg = SimulationConfig.from_dict(raw).with_overrides(
        {
            "time": {"duration": {"value": 2.0e-5, "unit": "s"}},
            "hydrodynamics": {"enabled": True},
            "output": {"archive_interval_s": 1.0e-5},
        }
    )
    simulator = Simulator(cfg)
    simulator.run(cfg.time.duration_s, output_policy="compact")
    assert len(simulator.hydro_samples) >= 2
    assert simulator.hydro_samples[0].t_s == 0.0
    assert (
        simulator.hydro_samples[0].positions_m.shape
        == simulator.hydro_samples[0].total_forces_N.shape
    )


def test_issue225_campaign_has_nine_fixed_attachment_conditions() -> None:
    campaign = load_yaml(
        Path("conf/phase2_multi_run/2010_project_hydrodynamics_issue225.yaml")
    )
    conditions = build_campaign_conditions(normalize_campaign_config(campaign))
    assert len(conditions) == 9
    assert campaign["base_overrides"]["seed.attach_seed"] == 0


def _state(t_s: float, positions_um: np.ndarray) -> SimulationState:
    return SimulationState(
        t=t_s,
        position_um=(0.0, 0.0, 0.0),
        quaternion=(0.0, 0.0, 0.0, 1.0),
        velocity_um_s=(1.0, 0.0, 0.0),
        omega_rad_s=(0.0, 0.0, 1.0),
        bead_positions_um=positions_um,
        flag_states=(),
        reverse_flagella=(),
    )


def test_campaign_analysis_writes_qc_filtered_force_flow_snapshots(
    tmp_path: Path,
) -> None:
    conditions = []
    positions_m = 1.0e-6 * np.asarray(
        [[-1.0, 1.0, 0.0], [-1.0, -1.0, 0.0], [1.0, 0.0, 0.0], [3.0, 0.0, 0.0]]
    )
    for n_flagella in (1, 2, 3):
        condition_id = f"n{n_flagella}__phase0"
        condition_dir = tmp_path / condition_id
        condition_dir.mkdir()
        positions_um = positions_m * 1e6
        save_state_archive(
            condition_dir / "state_archive.npz",
            [_state(0.0, positions_um), _state(0.001, positions_um)],
        )
        save_hydro_archive(
            condition_dir / "hydro_archive.npz",
            [
                HydroSample(
                    0.0,
                    positions_m,
                    np.asarray(
                        [
                            [1e-12, 0.0, 0.0],
                            [0.0, 1e-12, 0.0],
                            [0.0, -1e-12, 0.0],
                            [-1e-12, 0.0, 0.0],
                        ]
                    ),
                ),
                HydroSample(
                    0.001,
                    positions_m,
                    np.asarray(
                        [
                            [1e-12, 0.0, 0.0],
                            [0.0, 1e-12, 0.0],
                            [0.0, -1e-12, 0.0],
                            [-1e-12, 0.0, 0.0],
                        ]
                    ),
                ),
            ],
            bead_is_body=np.asarray([True, True, True, False]),
            bead_flagella_id=np.asarray([-1, -1, -1, 0]),
            bead_radius_m=1e-6,
            viscosity_Pa_s=1e-3,
            provenance={"hydrodynamics": {"model": "free_space_rpy"}},
        )
        nonbody_pass = n_flagella != 3
        body_fail = n_flagella == 3
        (condition_dir / "run_summary.json").write_text(
            '{"execution":{"status":"completed"},"gates":{"shape_nonbody":{"final_pass":'
            + str(nonbody_pass).lower()
            + '},"shape_body":{"any_fail":'
            + str(body_fail).lower()
            + "}}}"
        )
        conditions.append(
            {
                "condition_id": condition_id,
                "output_dir": str(condition_dir),
                "axis_values": {"n_flagella": n_flagella, "phase_seed": 0},
            }
        )
    (tmp_path / "run_manifest.json").write_text(
        __import__("json").dumps({"conditions": conditions})
    )
    output = analyze_campaign(tmp_path)
    manifest = __import__("json").loads((output / "analysis_manifest.json").read_text())
    assert len(manifest["included_conditions"]) == 2
    assert manifest["excluded_conditions"][0]["condition_id"] == "n3__phase0"
    assert (
        manifest["excluded_conditions"][0]["exclusion_reason"]
        == "strict_shape_qc_not_passed"
    )
    assert not (output / "hydrodynamics_comparison.csv").exists()
    assert not (output / "flagella_count_comparison.png").exists()
    assert not (output / "conditions").exists()
    assert manifest["layout"].startswith("3 phase-seed rows")


def test_body_fixed_flow_slice_uses_requested_dense_grid(tmp_path: Path) -> None:
    body = np.asarray(
        [
            [-1.0, 1.0, 0.0],
            [-1.0, -0.5, 0.8660254],
            [-1.0, -0.5, -0.8660254],
            [1.0, 1.0, 0.0],
            [1.0, -0.5, 0.8660254],
            [1.0, -0.5, -0.8660254],
        ]
    )
    path = tmp_path / "hydro_archive.npz"
    save_hydro_archive(
        path,
        [HydroSample(0.0, body * 1.0e-6, np.tile([1.0e-12, 0.0, 0.0], (len(body), 1)))],
        bead_is_body=np.ones(len(body), dtype=bool),
        bead_flagella_id=-np.ones(len(body), dtype=int),
        bead_radius_m=1.0e-6,
        viscosity_Pa_s=1.0e-3,
        provenance={},
    )
    archive = load_hydro_archive(path)
    points, velocity, frame = body_fixed_flow_slice(archive, 0)
    assert points.shape == (FLOW_SLICE_GRID_SIZE**2, 3)
    assert velocity.shape == points.shape
    np.testing.assert_allclose(frame, np.eye(3), atol=1.0e-7)


def test_hydrodynamics_replay_module_exposes_its_cli() -> None:
    result = subprocess.run(
        [sys.executable, "-m", "sim_swim.analysis.hydrodynamics_replay", "--help"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0
    assert "--row-axis" in result.stdout
