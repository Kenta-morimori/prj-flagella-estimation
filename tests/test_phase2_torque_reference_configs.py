from __future__ import annotations

import math
from pathlib import Path

import pytest

from sim_swim.analysis.multi_run_campaign import load_yaml
from sim_swim.sim.params import SimulationConfig


CONFIG_DIR = Path("conf/phase2_multi_run")


@pytest.mark.parametrize("config_path", sorted(CONFIG_DIR.glob("*.yaml")))
def test_standard_multi_run_profiles_keep_torque_and_reference_synced(
    config_path: Path,
) -> None:
    raw = load_yaml(config_path)
    if raw.get("kind") != "generic_multi_run":
        return
    base_cfg = SimulationConfig.from_dict(load_yaml(Path(raw["base_config"])))
    cfg = base_cfg.with_overrides(dict(raw.get("base_overrides") or {}))
    if cfg.is_motor_off_torque or cfg.motor.allow_reference_torque_mismatch:
        return
    assert math.isclose(
        abs(cfg.motor_torque_Nm),
        abs(cfg.reference_torque_Nm),
        rel_tol=1.0e-12,
        abs_tol=1.0e-30,
    ), config_path
