from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from sim_swim.analysis.flagella_count_behavior import save_state_archive
from sim_swim.analysis.motion_feature_study import (
    MotionFeatureStudyConfig,
    _axis_velocity,
    _basis,
    analyze_motion_feature_study,
)
from sim_swim.sim.core import SimulationState


def _csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _fixture(root: Path) -> None:
    conditions = []
    for n in (1, 4):
        sample = f"as000__ps000__nf{n:02d}"
        raw = root / sample
        raw.mkdir(parents=True)
        states = [
            SimulationState(
                t=index * 0.5,
                position_um=(index * 0.5, 0.0, 0.0),
                quaternion=(0.0, 0.0, 0.0, 1.0),
                velocity_um_s=(1.0, 0.0, 0.0),
                omega_rad_s=(0.0, 0.0, 0.0),
                bead_positions_um=np.asarray([[0.0, 0.0, 0.0]]),
            )
            for index in range(3)
        ]
        save_state_archive(raw / "state_archive.npz", states)
        _csv(
            raw / "step_summary.csv",
            [
                {
                    "t_s": state.t,
                    "finite_pass": True,
                    "shape_pass_nonbody_strict": True,
                    "first_fail_category_nonbody_strict": "none",
                }
                for state in states
            ],
        )
        axis_rows = []
        for state in states:
            for flag_id in range(n):
                axis_rows.append(
                    {
                        "t_s": state.t,
                        "flag_id": flag_id,
                        "axis_dir_x": -1.0,
                        "axis_dir_y": 0.0,
                        "axis_dir_z": 0.0,
                    }
                )
        _csv(raw / "flag_helix_axis_diagnostics.csv", axis_rows)
        conditions.append(
            {
                "condition_id": sample,
                "output_dir": str(raw),
                "axis_values": {"n_flagella": n, "attach_seed": 0, "phase_seed": 0},
            }
        )
    (root / "run_manifest.json").write_text(
        json.dumps({"conditions": conditions}), encoding="utf-8"
    )


@pytest.mark.light
def test_motion_feature_study_raw_campaign_writes_feature_contract(
    tmp_path: Path,
) -> None:
    run_dir, output_dir = tmp_path / "run", tmp_path / "analysis"
    _fixture(run_dir)
    analyze_motion_feature_study(
        MotionFeatureStudyConfig(
            run_dir=run_dir,
            output_dir=output_dir,
            durations_s=(1.0,),
            frame_rate_hz=2.0,
        )
    )
    rows = list(csv.DictReader((output_dir / "window_features_2d.csv").open()))
    assert rows and {row["n_flagella"] for row in rows} == {"1", "4"}
    assert all(float(row["mean_speed_um_s"]) == pytest.approx(1.0) for row in rows)
    assert (
        json.loads((output_dir / "manifest.json").read_text())["observability"][
            "2d_body_axis_angular_velocity"
        ]
        == "pixel_observable"
    )
    assert list((output_dir / "plots").glob("*.png"))


@pytest.mark.light
def test_motion_feature_axes_are_axial_and_projection_basis_is_orthonormal() -> None:
    vectors = np.asarray([[1.0, 0.0], [-1.0, 0.0]])
    assert _axis_velocity(vectors, np.asarray([0.0, 1.0]))[1] == pytest.approx(0.0)
    basis = _basis(((2.0, 0.0, 0.0), (1.0, 3.0, 0.0)))
    assert basis @ basis.T == pytest.approx(np.eye(2))
