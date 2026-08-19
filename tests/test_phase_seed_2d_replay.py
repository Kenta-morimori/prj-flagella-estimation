from __future__ import annotations

import json

import cv2
import numpy as np
import pytest

from sim_swim.analysis.flagella_count_behavior import save_state_archive
from sim_swim.analysis.phase_seed_2d_replay import (
    PhaseSeedReplayConfig,
    render_phase_seed_2d_replay,
)
from sim_swim.sim.core import SimulationState


def _states() -> list[SimulationState]:
    return [
        SimulationState(
            t=index * 0.5,
            position_um=(0.0, 0.0, 0.0),
            quaternion=(0.0, 0.0, 0.0, 1.0),
            velocity_um_s=(0.0, 0.0, 0.0),
            omega_rad_s=(0.0, 0.0, 0.0),
            bead_positions_um=np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
        )
        for index in range(3)
    ]


def _write_condition(root, *, n_flagella: int, phase_seed: int) -> None:
    condition_dir = root / f"as000__ps{phase_seed:03d}__nf{n_flagella:02d}"
    condition_dir.mkdir(parents=True)
    save_state_archive(condition_dir / "state_archive.npz", _states())
    first_fail = 0.25 if (n_flagella, phase_seed) == (4, 1) else None
    condition_dir.joinpath("run_summary.json").write_text(
        json.dumps(
            {
                "execution": {"status": "completed"},
                "gates": {
                    "shape_nonbody": {"first_observed_fail_t_s": first_fail},
                    "shape_body": {"first_observed_fail_t_s": first_fail},
                },
            }
        ),
        encoding="utf-8",
    )


@pytest.mark.light
def test_phase_seed_2d_replay_writes_fixed_grid_and_qc_panel(tmp_path) -> None:
    run_dir = tmp_path / "run"
    for n_flagella in range(1, 5):
        for phase_seed in range(3):
            _write_condition(run_dir, n_flagella=n_flagella, phase_seed=phase_seed)

    output_dir = render_phase_seed_2d_replay(
        PhaseSeedReplayConfig(
            run_dir=run_dir,
            output_dir=tmp_path / "replay",
            duration_s=1.0,
            fps=2.0,
            panel_width_px=120,
            panel_height_px=100,
        )
    )

    manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["frame_count"] == 2
    assert manifest["layout"] == {
        "rows": 4,
        "columns": 3,
        "row": "n_flagella",
        "column": "phase_seed",
    }
    assert len(manifest["conditions"]) == 12
    assert manifest["conditions"][-2]["condition_id"] == "as000__ps001__nf04"
    assert manifest["conditions"][-2]["first_fail_t_s"] == pytest.approx(0.25)
    assert (output_dir / "phase_seed_2d_grid.mp4").is_file()
    assert (output_dir / "run.log").is_file()

    image = cv2.imread(str(output_dir / "phase_seed_2d_grid_final.png"))
    assert image is not None
    # Row n=4, column ps=1 has entered QC failure by the final frame.
    assert tuple(image[3 * 100 + 1, 120 + 60]) == (225, 235, 255)
