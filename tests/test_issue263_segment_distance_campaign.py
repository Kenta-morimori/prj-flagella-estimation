"""Prepared nf6 distance controls preserve the archived compact model inputs."""

from copy import deepcopy
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pytest

from sim_swim.analysis.model_development_evaluation import _development_contract
from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    load_yaml,
    normalize_campaign_config,
)
from sim_swim.analysis.parallel_job import load_parallel_job, resolve_execution
from sim_swim.model.builder import ModelBuilder
from sim_swim.sim.params import SimulationConfig

pytestmark = pytest.mark.light
ROOT = Path(__file__).resolve().parents[1]
IDS = (
    "nf06__slots012345__distance_fixed_dt4us",
    "nf06__slots012345__distance_fixed_dt2us",
)


@pytest.mark.parametrize(
    ("tag", "stage", "steps"),
    [
        ("1tau", "short_screen", (10000, 20000)),
        ("25tau", "long_duration", (250000, 500000)),
    ],
)
def test_prepared_distance_campaign_preserves_legacy_geometry_and_physics(
    tag: str, stage: str, steps: tuple[int, int]
) -> None:
    path = (
        ROOT
        / f"conf/phase2_multi_run/2010_hex_nf6_segment_distance_{tag}_issue263.yaml"
    )
    campaign = normalize_campaign_config(load_yaml(path))
    contract = _development_contract(campaign)
    assert contract["stage"] == stage
    assert contract["expected_condition_count"] == 2
    assert contract["replay"]["camera_3d"] == "fixed"
    conditions = build_campaign_conditions(campaign)
    assert tuple(r["condition_id"] for r in conditions) == IDS
    old = normalize_campaign_config(
        load_yaml(
            ROOT
            / "conf/phase2_multi_run/2010_hex_compact_count13_body_contact_2s_issue263.yaml"
        )
    )
    legacy = next(
        r
        for r in build_campaign_conditions(old)
        if r["condition_id"] == "nf06__slots012345__bfoff"
    )
    base = SimulationConfig.from_dict(load_yaml(ROOT / campaign["base_config"]))
    reference = base.with_overrides(legacy["config_overrides"])
    reference_data = asdict(reference)
    time_fields = (
        "duration_s",
        "duration_unit",
        "duration_value",
        "dt_star",
        "integration_dt_star",
    )
    for field in time_fields:
        reference_data["time"].pop(field)
    for row, expected_steps, dt in zip(conditions, steps, (4e-6, 2e-6), strict=True):
        cfg = base.with_overrides(row["config_overrides"])
        assert cfg.dt_s == pytest.approx(dt)
        assert cfg.total_steps == expected_steps
        comparable = asdict(cfg)
        for field in time_fields:
            comparable["time"].pop(field)
        assert comparable == reference_data
        np.testing.assert_array_equal(
            ModelBuilder(cfg).build().positions_m,
            ModelBuilder(reference).build().positions_m,
        )
    job = load_parallel_job(
        ROOT
        / f"conf/phase2_parallel/issue263_initial_bundle/nf6_segment_distance_{tag}_job.yaml"
    )
    execution = resolve_execution(job, None)
    assert job.condition_ids == IDS
    assert execution.worker_policy == "cs10_qualified"
    assert execution.max_workers == 2
    assert set(execution.thread_environment.values()) == {"1"}
    assert campaign["metadata"]["historical_control"]["new_execution_required"] is False


def test_longer_distance_stage_changes_only_observation_duration() -> None:
    conditions = []
    for tag in ("1tau", "25tau"):
        campaign = normalize_campaign_config(
            load_yaml(
                ROOT
                / f"conf/phase2_multi_run/2010_hex_nf6_segment_distance_{tag}_issue263.yaml"
            )
        )
        conditions.append(build_campaign_conditions(campaign))
    for short, longer in zip(*conditions, strict=True):
        a, b = deepcopy(short["config_overrides"]), deepcopy(longer["config_overrides"])
        assert a["time"].pop("duration") == {"value": 1.0, "unit": "tau"}
        assert b["time"].pop("duration") == {"value": 25.0, "unit": "tau"}
        assert a == b
