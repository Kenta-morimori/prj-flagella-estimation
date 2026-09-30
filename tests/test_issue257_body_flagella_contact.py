from __future__ import annotations

from pathlib import Path

from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    load_yaml,
    normalize_campaign_config,
)


ROOT = Path(__file__).resolve().parents[1]
HEX_SCREEN = (
    ROOT
    / "conf/phase2_multi_run/2010_hex_project_body_flagella_contact_screen_issue257.yaml"
)
HEX_LONG = (
    ROOT
    / "conf/phase2_multi_run/2010_hex_project_body_flagella_contact_2s_issue257.yaml"
)
PROJECT_SCREEN = (
    ROOT
    / "conf/phase2_multi_run/2010_project_body_flagella_contact_screen_issue257.yaml"
)
PROJECT_LONG = (
    ROOT / "conf/phase2_multi_run/2010_project_body_flagella_contact_2s_issue257.yaml"
)


def _conditions(path: Path) -> tuple[dict, list[dict]]:
    config = normalize_campaign_config(load_yaml(path))
    return config, build_campaign_conditions(config)


def test_issue257_hex_contract_reuses_all_issue245_patterns_in_paired_arms() -> None:
    screen, screen_conditions = _conditions(HEX_SCREEN)
    long, long_conditions = _conditions(HEX_LONG)

    assert len(screen_conditions) == len(long_conditions) == 26
    assert screen["development_evaluation"]["expected_condition_count"] == 26
    assert long["development_evaluation"]["expected_condition_count"] == 26
    assert {
        condition["axis_values"]["body_flagella_repulsion"]
        for condition in screen_conditions
    } == {True, False}
    assert {
        condition["axis_values"]["attachment_pattern"]
        for condition in screen_conditions
    } == {
        "nf01__slots0",
        "nf02__slots01",
        "nf02__slots02",
        "nf02__slots03",
        "nf03__slots012",
        "nf03__slots013",
        "nf03__slots014",
        "nf03__slots024",
        "nf04__slots0123",
        "nf04__slots0124",
        "nf04__slots0134",
        "nf05__slots01234",
        "nf06__slots012345",
    }
    assert long["development_evaluation"]["stall_diagnostic"] == {
        "window_ms": 40,
        "relative_median_threshold": 0.1,
        "required_signals": ["body_speed_um_s", "body_roll_rate_hz"],
        "context_signal": "body_axis_step_angle_deg",
    }


def test_issue257_project_contract_is_supplemental_n4_to_n6_paired_screen() -> None:
    screen, screen_conditions = _conditions(PROJECT_SCREEN)
    long, long_conditions = _conditions(PROJECT_LONG)

    assert len(screen_conditions) == len(long_conditions) == 6
    assert {
        condition["axis_values"]["n_flagella"] for condition in screen_conditions
    } == {4, 5, 6}
    assert {
        condition["axis_values"]["body_flagella_repulsion"]
        for condition in screen_conditions
    } == {True, False}
    assert screen["metadata"]["canonical"] is False
    assert long["metadata"]["stall_diagnostic"]["window_ms"] == 40
