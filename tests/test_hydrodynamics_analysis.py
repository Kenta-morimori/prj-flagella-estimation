from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from sim_swim.analysis.hydrodynamics import (
    HYDRO_ARCHIVE_FORMAT,
    HydroSample,
    load_hydro_archive,
    rpy_flow_velocity,
    save_hydro_archive,
    velocity_contributions,
)
from sim_swim.analysis.hydrodynamics_campaign import (
    _body_angular_speed,
    _body_frame,
    _comparison_series,
    analyze_campaign,
)
from sim_swim.analysis.hydrodynamics_replay import (
    campaign_nflagella_phase0_condition_ids,
    overlay_manifest,
)
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


def test_campaign_analysis_writes_comparison_and_full_run_slice(tmp_path: Path) -> None:
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
        (condition_dir / "run_summary.json").write_text(
            '{"execution":{"status":"completed"},"gates":{"shape_nonbody":{"final_pass":true},"shape_body":{"any_fail":false}}}'
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
    assert (output / "hydrodynamics_comparison.csv").is_file()
    assert (output / "flagella_count_comparison.png").is_file()
    assert (output / "body_fixed_axial_flow.png").is_file()


def test_body_frame_and_angular_speed_include_axial_roll() -> None:
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
    quarter_roll = np.asarray([[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]])
    frames = [
        _body_frame(body, np.ones(len(body), dtype=bool)),
        _body_frame(body @ quarter_roll.T, np.ones(len(body), dtype=bool)),
    ]
    np.testing.assert_allclose(frames[0], np.eye(3), atol=1.0e-7)
    assert _body_angular_speed(frames, np.asarray([0.0, 0.5]), 1) == pytest.approx(
        np.pi, rel=1.0e-7
    )


def test_campaign_overlay_contract_is_fixed_camera_three_panels() -> None:
    payload = overlay_manifest(["n1__phase0", "n2__phase0", "n3__phase0"], fps=25.0)
    assert payload["follow_camera_3d"] is False
    assert payload["grid_shape"] == [3, 3, 3]
    assert payload["panel_count"] == 3


def test_hydrodynamics_replay_module_exposes_its_cli() -> None:
    result = subprocess.run(
        [sys.executable, "-m", "sim_swim.analysis.hydrodynamics_replay", "--help"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0
    assert "--campaign-nflagella-phase0" in result.stdout


def test_campaign_overlay_selects_issue215_style_condition_ids(tmp_path: Path) -> None:
    conditions = [
        {
            "condition_id": f"as000__ps000__nf{count:02d}",
            "axis_values": {"attach_seed": 0, "phase_seed": 0, "n_flagella": count},
        }
        for count in (1, 2, 3)
    ] + [
        {
            "condition_id": "as001__ps000__nf01",
            "axis_values": {"attach_seed": 1, "phase_seed": 0, "n_flagella": 1},
        }
    ]
    (tmp_path / "run_manifest.json").write_text(
        __import__("json").dumps({"conditions": conditions})
    )
    assert campaign_nflagella_phase0_condition_ids(tmp_path) == [
        "as000__ps000__nf01",
        "as000__ps000__nf02",
        "as000__ps000__nf03",
    ]


def test_comparison_series_separates_attachment_seeds() -> None:
    rows = [
        {"attach_seed": attach, "phase_seed": 0, "n_flagella": count}
        for attach in (0, 1)
        for count in (1, 2)
    ]
    series = _comparison_series(rows)
    assert [label for label, _ in series] == ["attach 0, phase 0", "attach 1, phase 0"]
    assert [[row["n_flagella"] for row in values] for _, values in series] == [
        [1, 2],
        [1, 2],
    ]
