from __future__ import annotations

import json
from pathlib import Path

import pytest

from sim_swim.analysis.flagella_count_behavior import validate_replay_fps
from sim_swim.sim.core import SimulationInterrupted, Simulator
from sim_swim.sim.params import SimulationConfig


def _cfg(*, policy: str, checkpoint_interval_steps: int = 2500) -> SimulationConfig:
    raw = {
        "time": {"duration_s": 0.002, "dt_s": 0.001, "dt_star": 0.001},
        "motor": {"torque_Nm": 1.0e-21, "reference_torque_Nm": 1.0e-21},
        "output": {
            "policy": policy,
            "archive_interval_s": 0.001,
            "checkpoint_interval_steps": checkpoint_interval_steps,
        },
        "flagella": {"n_flagella": 0},
    }
    return SimulationConfig.from_dict(raw)


def test_compact_keeps_every_step_qc_without_step_csv(tmp_path: Path) -> None:
    cfg = _cfg(policy="compact")
    assert cfg.output.checkpoint_interval_steps == 2500
    states = Simulator(cfg).run(cfg.time.duration_s, step_summary_dir=tmp_path)
    summary = json.loads((tmp_path / "run_summary.json").read_text())
    assert not (tmp_path / "step_summary.csv").exists()
    assert summary["execution"]["row_count"] == cfg.total_steps
    assert summary["gates"]["finite"]["status"] == "available"
    assert summary["gates"]["shape_body"]["status"] == "available"
    assert states[0].t == 0.0 and states[-1].t == pytest.approx(cfg.time.duration_s)
    performance = json.loads((tmp_path / "performance.json").read_text())
    assert performance["saved_state_count"] == len(states)
    assert performance["steps_per_s"] > 0.0


def test_compact_archive_fps_limit_and_debug_csv_compatibility(tmp_path: Path) -> None:
    compact_cfg = _cfg(policy="compact")
    states = Simulator(compact_cfg).run(
        compact_cfg.time.duration_s, step_summary_dir=tmp_path / "compact"
    )
    validate_replay_fps(states, 1000.0)
    with pytest.raises(ValueError, match="exceeds archive density"):
        validate_replay_fps(states, 1001.0)
    debug_cfg = _cfg(policy="debug")
    Simulator(debug_cfg).run(
        debug_cfg.time.duration_s, step_summary_dir=tmp_path / "debug"
    )
    assert (tmp_path / "debug" / "step_summary.csv").is_file()
    debug_summary = json.loads((tmp_path / "debug" / "run_summary.json").read_text())
    assert debug_summary["execution"] == {
        "status": "completed",
        "row_count": debug_cfg.total_steps,
        "step_indices_contiguous_from_zero": True,
        "observed_first_t_s": pytest.approx(debug_cfg.dt_s),
        "observed_final_t_s": pytest.approx(debug_cfg.final_state_t_s),
        "expected_total_steps": debug_cfg.total_steps,
        "expected_final_step_summary_t_s": pytest.approx(debug_cfg.final_state_t_s),
        "reason": "matched time manifest",
    }


def test_compact_archive_interpolates_to_the_requested_time_grid() -> None:
    cfg = SimulationConfig.from_dict(
        {
            "time": {"duration_s": 0.003, "dt_s": 0.0001, "dt_star": 0.0001},
            "motor": {"torque_Nm": 1.0e-21, "reference_torque_Nm": 1.0e-21},
            "output": {"policy": "compact", "archive_interval_s": 0.00015},
            "flagella": {"n_flagella": 0},
        }
    )

    states = Simulator(cfg).run(cfg.time.duration_s)

    assert [state.t for state in states] == pytest.approx(
        [0.00015 * index for index in range(21)]
    )


def test_compact_exception_writes_partial_summary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    cfg = _cfg(policy="compact")
    sim = Simulator(cfg)
    original = sim.engine.step
    calls = 0

    def fail_after_first(dt_star: float):
        nonlocal calls
        calls += 1
        if calls > 1:
            raise RuntimeError("intentional")
        return original(dt_star)

    monkeypatch.setattr(sim.engine, "step", fail_after_first)
    with pytest.raises(RuntimeError, match="intentional"):
        sim.run(cfg.time.duration_s, step_summary_dir=tmp_path)
    summary = json.loads((tmp_path / "run_summary.json").read_text())
    assert summary["execution"]["status"] == "partial"
    assert "RuntimeError" in summary["execution"]["reason"]


def test_compact_checkpoint_callback_reports_raw_boundary_evidence(
    tmp_path: Path,
) -> None:
    cfg = _cfg(policy="compact", checkpoint_interval_steps=1)
    checkpoints: list[dict[str, object]] = []

    Simulator(cfg).run(
        cfg.time.duration_s,
        step_summary_dir=tmp_path,
        checkpoint_callback=lambda **payload: checkpoints.append(payload),
    )

    assert [checkpoint["completed_steps"] for checkpoint in checkpoints] == [1, 2, 2]
    assert [checkpoint["status"] for checkpoint in checkpoints] == [
        "running",
        "running",
        "completed",
    ]
    assert all(checkpoint["diagnostic_row"] is not None for checkpoint in checkpoints)
    assert all(checkpoint["body_row"] is not None for checkpoint in checkpoints)


def test_compact_cooperative_interrupt_writes_partial_summary_and_checkpoint(
    tmp_path: Path,
) -> None:
    cfg = _cfg(policy="compact", checkpoint_interval_steps=2)
    checkpoints: list[dict[str, object]] = []

    with pytest.raises(SimulationInterrupted, match="test interruption"):
        Simulator(cfg).run(
            cfg.time.duration_s,
            step_summary_dir=tmp_path,
            checkpoint_callback=lambda **payload: checkpoints.append(payload),
            interrupt_requested=lambda: "test interruption",
        )

    assert checkpoints[-1]["status"] == "partial"
    assert checkpoints[-1]["completed_steps"] == 1
    summary = json.loads((tmp_path / "run_summary.json").read_text())
    performance = json.loads((tmp_path / "performance.json").read_text())
    assert summary["execution"]["status"] == "partial"
    assert performance["completed_steps"] == 1
