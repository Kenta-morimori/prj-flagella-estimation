from __future__ import annotations

import json
from pathlib import Path

import pytest

from sim_swim.analysis.flagella_count_behavior import validate_replay_fps
from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig


def _cfg(*, policy: str) -> SimulationConfig:
    raw = {
        "time": {"duration_s": 0.002, "dt_s": 0.001, "dt_star": 0.001},
        "output": {"policy": policy, "archive_interval_s": 0.001},
        "flagella": {"n_flagella": 0},
    }
    return SimulationConfig.from_dict(raw)


def test_compact_keeps_every_step_qc_without_step_csv(tmp_path: Path) -> None:
    cfg = _cfg(policy="compact")
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
