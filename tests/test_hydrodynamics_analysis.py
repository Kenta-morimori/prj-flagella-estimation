from pathlib import Path

import numpy as np

from sim_swim.analysis.hydrodynamics import (
    HYDRO_ARCHIVE_FORMAT,
    HydroSample,
    load_hydro_archive,
    rpy_flow_velocity,
    save_hydro_archive,
    velocity_contributions,
)
from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    load_yaml,
    normalize_campaign_config,
)
from sim_swim.sim.core import Simulator
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
