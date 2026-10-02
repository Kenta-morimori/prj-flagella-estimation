from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from sim_swim.analysis.model_development_evaluation import (
    _without_override_paths,
    _write_replay_input,
    collect_rows,
)
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
    assert long["development_evaluation"]["reused_source_campaigns"] == {
        "conf/phase2_multi_run/2010_hex_project_long_duration_2s_issue245.yaml": {
            "axis_defaults": {"body_flagella_repulsion": True},
            "relaxed_config_override_paths": [
                "potentials.spring_spring_repulsion.body_flagella_enabled"
            ],
        }
    }


def _summary() -> dict:
    metrics = {
        "local_attach_first_rel_err": 0.01,
        "hook_len_rel_err_max": 0.01,
        "flag_bond_rel_err_max": 0.01,
        "flag_bend_err_max_deg": 1.0,
        "flag_torsion_err_max_deg": 1.0,
        "body_spring_max_stretch_ratio": 0.01,
        "motor_force_balance_residual_ratio": 1e-16,
        "motor_torque_balance_residual_ratio": 0.01,
    }
    return {
        "execution": {"status": "completed"},
        "gates": {
            "finite": {"any_fail": False},
            "shape_body": {"any_fail": False},
            "shape_nonbody": {"any_fail": False},
        },
        "all_step_metrics": {name: {"max": value} for name, value in metrics.items()},
    }


def _write_source(
    root: Path,
    *,
    config: dict,
    conditions: list[dict],
    campaign_config: str,
    legacy_on: bool,
) -> Path:
    root.mkdir(parents=True)
    records: list[dict] = []
    rows: list[dict[str, str]] = []
    for expected in conditions:
        source_id = str(expected["condition_id"])
        record = json.loads(json.dumps(expected))
        if legacy_on:
            source_id = source_id.removesuffix("__bfon")
            record["condition_id"] = source_id
            record["axis_values"].pop("body_flagella_repulsion")
            record["axis_labels"].pop("body_flagella_repulsion")
            record["config_overrides"]["potentials"]["spring_spring_repulsion"].pop(
                "body_flagella_enabled"
            )
        output = root / "conditions" / source_id
        output.mkdir(parents=True)
        (output / "run_summary.json").write_text(json.dumps(_summary()))
        (output / "performance.json").write_text(
            json.dumps({"wall_time_s": 1.0, "steps_per_s": 10.0})
        )
        (output / "state_archive.npz").write_bytes(b"complete archive")
        (output / "initial_geometry_summary.json").write_text(
            json.dumps(
                {
                    "geometry": {
                        "body_slots_per_layer": 6,
                        "body_layers": 5,
                        "body_beads": 30,
                        "attachment_topology": [
                            {
                                "flag_id": index,
                                "body_bead_index": 12 + slot,
                                "layer": 2,
                                "slot": slot,
                            }
                            for index, slot in enumerate(
                                record["axis_values"]["attachment_slots"]
                            )
                        ],
                    }
                }
            )
        )
        record["output_dir"] = str(output)
        records.append(record)
        rows.append({"condition_id": source_id, "completed": "True"})
    with (root / "summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["condition_id", "completed"])
        writer.writeheader()
        writer.writerows(rows)
    (root / "run_manifest.json").write_text(
        json.dumps(
            {
                "base_config": "conf/sim_swim_2010_hex.yaml",
                "campaign_config": campaign_config,
                "model_profile": load_yaml(ROOT / config["base_config"])[
                    "model_profile"
                ],
                "git": {"commit": "f41a693", "is_clean": True},
                "conditions": records,
            }
        )
    )
    return root


def test_issue257_long_evaluation_reuses_issue245_on_with_off_source(
    tmp_path: Path,
) -> None:
    config = load_yaml(HEX_LONG)
    _, conditions = _conditions(HEX_LONG)
    on = _write_source(
        tmp_path / "issue245_on",
        config=config,
        conditions=[
            item for item in conditions if item["condition_id"].endswith("bfon")
        ],
        campaign_config="conf/phase2_multi_run/2010_hex_project_long_duration_2s_issue245.yaml",
        legacy_on=True,
    )
    off = _write_source(
        tmp_path / "issue257_off",
        config=config,
        conditions=[
            item for item in conditions if item["condition_id"].endswith("bfoff")
        ],
        campaign_config="conf/phase2_multi_run/2010_hex_project_body_flagella_contact_2s_issue257.yaml",
        legacy_on=False,
    )

    rows, provenance = collect_rows(config=config, run_dirs=[on, off])

    assert len(rows) == 26
    assert {row["condition_id"].split("__")[-1] for row in rows} == {
        "bfon",
        "bfoff",
    }
    assert sum(bool(row["source_reused"]) for row in rows) == 13
    assert {row["source_git_commit"] for row in rows} == {"f41a693"}
    assert len(provenance) == 2

    manifest = json.loads((on / "run_manifest.json").read_text())
    manifest["conditions"][0]["artifact_sha256"] = {"state_archive.npz": "wrong"}
    (on / "run_manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        collect_rows(config=config, run_dirs=[on, off])


def test_issue257_parallel_campaign_uses_parent_job_git_provenance(
    tmp_path: Path,
) -> None:
    config = load_yaml(HEX_LONG)
    _, conditions = _conditions(HEX_LONG)
    on = _write_source(
        tmp_path / "issue245_on",
        config=config,
        conditions=[
            item for item in conditions if item["condition_id"].endswith("bfon")
        ],
        campaign_config="conf/phase2_multi_run/2010_hex_project_long_duration_2s_issue245.yaml",
        legacy_on=True,
    )
    source = _write_source(
        tmp_path / "parallel" / "campaign",
        config=config,
        conditions=[
            item for item in conditions if item["condition_id"].endswith("bfoff")
        ],
        campaign_config="conf/phase2_multi_run/2010_hex_project_body_flagella_contact_2s_issue257.yaml",
        legacy_on=False,
    )
    manifest = json.loads((source / "run_manifest.json").read_text())
    manifest.pop("git")
    (source / "run_manifest.json").write_text(json.dumps(manifest))
    (source / "job_manifest.json").write_text(
        json.dumps(
            {
                "provenance": {
                    "git": {"commit": "e09c632", "status": "## HEAD (no branch)"}
                }
            }
        )
    )

    rows, provenance = collect_rows(config=config, run_dirs=[on, source])

    assert len(rows) == 26
    assert {row["source_git_commit"] for row in rows} == {"f41a693", "e09c632"}
    assert provenance[1]["git"]["source"] == "parallel_job_manifest"


def test_issue257_reused_selector_comparison_accepts_omitted_parent_mapping() -> None:
    assert (
        _without_override_paths(
            {
                "potentials": {
                    "spring_spring_repulsion": {"body_flagella_enabled": True}
                }
            },
            ["potentials.spring_spring_repulsion.body_flagella_enabled"],
        )
        == {}
    )


def test_issue257_replay_input_falls_back_from_stale_source_config_path(
    tmp_path: Path,
) -> None:
    config = load_yaml(HEX_LONG)
    _, conditions = _conditions(HEX_LONG)
    source = _write_source(
        tmp_path / "source",
        config=config,
        conditions=[
            item for item in conditions if item["condition_id"].endswith("bfoff")
        ],
        campaign_config="conf/phase2_multi_run/2010_hex_project_body_flagella_contact_2s_issue257.yaml",
        legacy_on=False,
    )
    manifest = json.loads((source / "run_manifest.json").read_text())
    manifest["source_config_path"] = "/missing/cs10/conf/sim_swim_2010_hex.yaml"
    (source / "run_manifest.json").write_text(json.dumps(manifest))

    replay_input = _write_replay_input(
        output_dir=tmp_path / "evaluation", run_dirs=[source], config=config
    )

    replay_manifest = json.loads((replay_input / "run_manifest.json").read_text())
    assert replay_manifest["base_config"] == "conf/sim_swim_2010_hex.yaml"


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
