from __future__ import annotations

from pathlib import Path

from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    load_yaml,
    normalize_campaign_config,
)
from sim_swim.sim.params import SimulationConfig


def test_hook_diagnostic_changes_only_local_hook_scale() -> None:
    config = normalize_campaign_config(
        load_yaml(
            Path(
                "conf/phase2_multi_run/"
                "reference_torque_2010_tracking_scale2_hook_diagnostic.yaml"
            )
        )
    )
    conditions = build_campaign_conditions(config)

    assert [row["axis_values"]["hook_scale"] for row in conditions] == [
        1.0,
        1.25,
        1.5,
        2.0,
    ]
    for condition in conditions:
        cfg = SimulationConfig.from_dict(
            load_yaml(Path(config["base_config"]))
        ).with_overrides(condition["config_overrides"])
        assert cfg.motor_torque_Nm == 5.0e-20
        assert cfg.reference_torque_Nm == 5.0e-20
        assert cfg.torque_for_forces_Nm == 5.0e-20
        assert cfg.time.duration_s == 0.02
        assert cfg.dt_star == 1.0e-4
