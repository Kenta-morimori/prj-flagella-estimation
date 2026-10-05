from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from sim_swim.analysis.model_development_evaluation import (
    SCREEN_METRICS,
    _plot_reaction_pairs,
)
from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    geometry_preflight,
    load_yaml,
    normalize_campaign_config,
)
from sim_swim.analysis.parallel_job import (
    build_plan,
    load_parallel_job,
    resolve_execution,
)
from sim_swim.model.builder import ModelBuilder
from sim_swim.sim.params import SimulationConfig

ROOT = Path(__file__).resolve().parents[1]
CAMPAIGNS = {
    "hex": ROOT
    / "conf/phase2_multi_run/2010_hex_project_motor_reaction_1tau_issue255.yaml",
    "project": ROOT
    / "conf/phase2_multi_run/2010_project_motor_reaction_1tau_issue255.yaml",
}


@pytest.mark.parametrize("model_name,expected", [("hex", 26), ("project", 6)])
def test_issue255_geometry_and_paired_contract(
    model_name: str, expected: int, tmp_path: Path
) -> None:
    campaign = normalize_campaign_config(load_yaml(CAMPAIGNS[model_name]))
    conditions = build_campaign_conditions(campaign)
    preflight = geometry_preflight(campaign, conditions)
    job = load_parallel_job(
        ROOT
        / f"conf/phase2_parallel/issue255_motor_reaction/{model_name}_1tau_job.yaml"
    )
    execution = resolve_execution(job, None)
    plan = build_plan(job, execution, tmp_path / "plan")
    assert (
        len(conditions)
        == len(preflight)
        == job.task_count
        == len(plan["configs"])
        == expected
    )
    assert job.preflight == "geometry_all_conditions"
    assert execution.max_workers == min(8, expected)
    assert execution.worker_policy == "cs10_qualified"
    assert set(execution.thread_environment.values()) == {"1"}
    assert len({item["output_dir"] for item in plan["configs"]}) == expected
    assert max(x["initial_hook_body_axis_error_deg"] for x in preflight.values()) < 1e-6
    assert (
        max(x["initial_hook_angle_max_abs_error_deg"] for x in preflight.values())
        < 1e-6
    )
    assert max(x["initial_hook_force_norm_N"] for x in preflight.values()) < 1e-18

    base = load_yaml(ROOT / campaign["base_config"])
    by_shape: dict[str, list[tuple[SimulationConfig, object]]] = {}
    for condition in conditions:
        assert isinstance(condition["axis_values"]["n_flagella"], int)
        assert float(condition["axis_values"]["motor_torque"]) == 2.5e-20
        cfg = SimulationConfig.from_dict(base).with_overrides(
            condition["config_overrides"]
        )
        assert cfg.flagella.initial_hook_force_neutral
        assert cfg.flagella.initial_hook_body_axis_perpendicular
        assert not cfg.potentials.spring_spring_repulsion.body_flagella_enabled
        assert not cfg.brownian.enabled and not cfg.motor.enable_switching
        assert cfg.time.duration_unit == "tau" and cfg.time.duration_value == 1.0
        assert cfg.time.integration_dt_star == 1e-4
        model = ModelBuilder(cfg).build()
        shape = condition["condition_id"].rsplit("__rx", 1)[0]
        by_shape.setdefault(shape, []).append((cfg, model))
        assert (
            preflight[condition["condition_id"]][
                "initial_min_nonattached_bead_distance_m"
            ]
            >= 2 * model.bead_radius_m
        )
        assert (
            preflight[condition["condition_id"]][
                "initial_min_attachment_bead_distance_m"
            ]
            >= 2 * model.bead_radius_m
        )

    assert len(by_shape) == expected // 2
    if model_name == "project":
        assert {pair[0][0].flagella.n_flagella for pair in by_shape.values()} == {
            1,
            2,
            3,
        }
    for pair in by_shape.values():
        assert len(pair) == 2
        (axis_cfg, axis_model), (full_cfg, full_model) = pair
        assert not axis_cfg.motor.body_reaction_full_vector
        assert full_cfg.motor.body_reaction_full_vector
        np.testing.assert_array_equal(axis_model.positions_m, full_model.positions_m)
        np.testing.assert_array_equal(
            axis_model.flagella_attach_body_indices,
            full_model.flagella_attach_body_indices,
        )

        old_cfg = axis_cfg.with_overrides(
            {
                "flagella": {
                    "initial_hook_force_neutral": False,
                    "initial_hook_body_axis_perpendicular": False,
                }
            }
        )
        old = ModelBuilder(old_cfg).build()
        np.testing.assert_array_equal(
            old.positions_m[old.body_indices],
            axis_model.positions_m[axis_model.body_indices],
        )
        for new_ids, old_ids in zip(axis_model.flagella_indices, old.flagella_indices):
            new_points = axis_model.positions_m[new_ids]
            old_points = old.positions_m[old_ids]
            np.testing.assert_allclose(new_points[0], old_points[0], atol=1e-20)
            np.testing.assert_allclose(
                np.linalg.norm(new_points[:, None] - new_points[None, :], axis=-1),
                np.linalg.norm(old_points[:, None] - old_points[None, :], axis=-1),
                atol=1e-15,
            )


def test_perpendicular_mode_requires_neutral_mode() -> None:
    base = load_yaml(ROOT / "conf/sim_swim_2010_hex.yaml")
    cfg = SimulationConfig.from_dict(base).with_overrides(
        {"flagella": {"initial_hook_body_axis_perpendicular": True}}
    )
    with pytest.raises(ValueError, match="requires.*initial_hook_force_neutral"):
        ModelBuilder(cfg).build()


def test_reaction_pair_heatmap_retains_both_arms(tmp_path: Path) -> None:
    rows = []
    for arm, status, residual in ((False, "fail", 0.2), (True, "pass", 1e-15)):
        row = {
            "n_flagella": 1,
            "attachment_pattern": "nf01__slots0",
            "body_reaction_full_vector": arm,
            "screen_status": status,
            **{metric: 0.0 for metric, _ in SCREEN_METRICS},
        }
        row["motor_torque_balance_residual_ratio"] = residual
        rows.append(row)
    output = tmp_path / "reaction_pairs.png"
    _plot_reaction_pairs(rows, output)
    assert output.is_file() and output.stat().st_size > 1000
    with pytest.raises(ValueError, match="exactly two arms"):
        _plot_reaction_pairs(rows[:1], output)


@pytest.mark.parametrize("model_name,expected_shapes", [("hex", 13), ("project", 3)])
def test_preview_adds_axial_projection_and_overview(
    model_name: str, expected_shapes: int, tmp_path: Path
) -> None:
    output_dir = tmp_path / model_name
    environment = {**os.environ, "MPLCONFIGDIR": str(tmp_path / "mpl_cache")}
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/01_simulate_swimming/preview_initial_geometry.py"),
            "--config",
            str(CAMPAIGNS[model_name]),
            "--output-dir",
            str(output_dir),
        ],
        cwd=ROOT,
        env=environment,
        check=True,
        capture_output=True,
        text=True,
    )
    manifest = json.loads((output_dir / "manifest.json").read_text())
    campaign = manifest["campaigns"][0]
    assert campaign["unique_shape_count"] == expected_shapes
    assert len(campaign["overview_condition_ids"]) == min(expected_shapes, 6)
    assert "onto y-z" in campaign["axial_projection"]
    for path_key, hash_key in (
        ("image", "image_sha256"),
        ("axial_image", "axial_image_sha256"),
        ("overview_image", "overview_image_sha256"),
    ):
        path = Path(campaign[path_key])
        assert path.is_file()
        assert hashlib.sha256(path.read_bytes()).hexdigest() == campaign[hash_key]
