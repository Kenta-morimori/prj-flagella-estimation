"""Contracts retained for executed #257 campaigns and the pending hook screen."""

from __future__ import annotations

from pathlib import Path

from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    load_yaml,
    normalize_campaign_config,
)
from sim_swim.analysis.parallel_job import load_parallel_job, resolve_execution

ROOT = Path(__file__).resolve().parents[1]
JOB_DIR = ROOT / "conf/phase2_parallel/issue257_body_flagella_contact"
CONF_DIR = ROOT / "conf/phase2_multi_run"


def _conditions(name: str) -> list[dict]:
    config = normalize_campaign_config(load_yaml(CONF_DIR / name))
    return build_campaign_conditions(config)


def test_executed_hex_campaigns_keep_13_patterns_and_off_only_long_run() -> None:
    screen = _conditions("2010_hex_project_body_flagella_contact_screen_issue257.yaml")
    long = _conditions("2010_hex_project_body_flagella_contact_2s_issue257.yaml")
    assert len(screen) == 26
    assert len(long) == 26
    assert {condition["axis_values"]["attachment_pattern"] for condition in screen} == {
        condition["axis_values"]["attachment_pattern"] for condition in long
    }
    assert {
        condition["axis_values"]["body_flagella_repulsion"] for condition in screen
    } == {True, False}
    off = [
        condition
        for condition in long
        if condition["axis_values"]["body_flagella_repulsion"] is False
    ]
    assert len(off) == 13
    job = load_parallel_job(JOB_DIR / "hex_2s_job.yaml")
    assert set(job.condition_ids) == {condition["condition_id"] for condition in off}
    assert job.preflight == "geometry_all_conditions"
    assert resolve_execution(job, None).worker_policy == "cs10_qualified"
    assert resolve_execution(job, None).max_workers == 8


def test_executed_project_screen_is_supplemental_only() -> None:
    screen = _conditions("2010_project_body_flagella_contact_screen_issue257.yaml")
    assert len(screen) == 6
    assert {condition["axis_values"]["n_flagella"] for condition in screen} == {
        4,
        5,
        6,
    }
    assert {
        condition["axis_values"]["body_flagella_repulsion"] for condition in screen
    } == {True, False}
    job = load_parallel_job(JOB_DIR / "project_screen_job.yaml")
    assert set(job.condition_ids) == {condition["condition_id"] for condition in screen}
