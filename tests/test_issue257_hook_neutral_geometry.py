from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

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
from sim_swim.dynamics.forces import compute_hook_forces
from sim_swim.model.builder import (
    ModelBuilder,
    _neutralize_initial_hooks_without_bead_overlap,
)
from sim_swim.sim.debug_summary import _triplet_angles_rad
from sim_swim.sim.params import SimulationConfig


ROOT = Path(__file__).resolve().parents[1]
CAMPAIGN = (
    ROOT / "conf/phase2_multi_run/2010_hex_project_hook_neutral_screen_issue257.yaml"
)
JOB = (
    ROOT
    / "conf/phase2_parallel/issue257_body_flagella_contact/hex_hook_neutral_screen_job.yaml"
)


@pytest.fixture(scope="module")
def screen():
    campaign = normalize_campaign_config(load_yaml(CAMPAIGN))
    return campaign, build_campaign_conditions(campaign)


def test_new_short_screen_is_opt_in_and_cs10_parallel(screen, tmp_path: Path) -> None:
    campaign, conditions = screen
    job = load_parallel_job(JOB)
    plan = build_plan(job, resolve_execution(job, None), tmp_path / "plan")
    assert len(conditions) == job.task_count == len(plan["configs"]) == 26
    assert set(job.condition_ids) == {item["condition_id"] for item in conditions}
    assert job.preflight == "geometry_all_conditions"
    assert resolve_execution(job, None).max_workers == 8
    assert resolve_execution(job, None).worker_policy == "cs10_qualified"
    assert len({item["output_dir"] for item in plan["configs"]}) == 26
    assert all(
        item["geometry_preflight"]["initial_hook_angle_max_abs_error_deg"] < 1e-6
        for item in plan["configs"]
    )
    assert all(
        item["geometry_preflight"]["initial_min_nonattached_bead_distance_m"] > 0
        for item in plan["configs"]
    )
    assert campaign["base_overrides"]["flagella"]["initial_hook_force_neutral"] is True


@pytest.mark.parametrize("phase_seed", [0, 1, 7])
def test_neutral_geometry_preserves_flagella_and_clears_beads(
    screen, phase_seed: int
) -> None:
    campaign, conditions = screen
    base = load_yaml(ROOT / campaign["base_config"])
    assert SimulationConfig.from_dict(base).flagella.initial_hook_force_neutral is False
    for condition in conditions:
        overrides = {**condition["config_overrides"], "seed.phase_seed": phase_seed}
        neutral_cfg = SimulationConfig.from_dict(base).with_overrides(overrides)
        old_cfg = SimulationConfig.from_dict(base).with_overrides(
            {**overrides, "flagella.initial_hook_force_neutral": False}
        )
        neutral = ModelBuilder(neutral_cfg).build()
        old = ModelBuilder(old_cfg).build()
        np.testing.assert_array_equal(
            neutral.positions_m[neutral.body_indices], old.positions_m[old.body_indices]
        )
        np.testing.assert_array_equal(
            neutral.flagella_attach_body_indices, old.flagella_attach_body_indices
        )
        np.testing.assert_array_equal(
            neutral.flagella_initial_phases_rad, old.flagella_initial_phases_rad
        )
        angles = np.degrees(
            _triplet_angles_rad(neutral.positions_m, neutral.hook_triplets)
        )
        np.testing.assert_allclose(angles, 90.0, atol=1e-6)
        force = compute_hook_forces(
            neutral.positions_m,
            neutral.hook_triplets,
            neutral_cfg.hook.kb_over_T * abs(neutral_cfg.motor.torque_Nm),
            neutral_cfg.hook.threshold_deg,
        )
        assert np.linalg.norm(force) < 1e-18
        body = neutral.positions_m[neutral.body_indices]
        body_center = body.mean(axis=0)
        flags = []
        for idx, attach_idx in zip(
            neutral.flagella_indices, neutral.flagella_attach_body_indices
        ):
            new_points = neutral.positions_m[idx]
            old_points = old.positions_m[idx]
            np.testing.assert_allclose(
                np.diff(new_points, axis=0), np.diff(old_points, axis=0), atol=1e-20
            )
            np.testing.assert_allclose(
                np.linalg.norm(new_points[0] - body[int(attach_idx)]),
                neutral_cfg.hook.length_over_b * neutral_cfg.b_m,
                atol=1e-15,
            )
            assert (
                np.dot(
                    new_points[0] - body[int(attach_idx)],
                    body[int(attach_idx)] - body_center,
                )
                > 0
            )
            assert (
                np.min(
                    np.linalg.norm(
                        new_points[:, None]
                        - np.delete(body, int(attach_idx), axis=0)[None, :],
                        axis=-1,
                    )
                )
                >= 2 * neutral.bead_radius_m - 1e-15
            )
            flags.append(new_points)
        for i in range(len(flags)):
            for j in range(i + 1, len(flags)):
                assert (
                    np.min(
                        np.linalg.norm(flags[i][:, None] - flags[j][None, :], axis=-1)
                    )
                    >= 2 * neutral.bead_radius_m - 1e-15
                )


def test_neutral_preflight_records_initial_geometry(screen) -> None:
    campaign, conditions = screen
    records = geometry_preflight(campaign, conditions)
    assert len(records) == 26
    assert (
        max(item["initial_hook_angle_max_abs_error_deg"] for item in records.values())
        < 1e-6
    )
    assert max(item["initial_hook_force_norm_N"] for item in records.values()) < 1e-18
    assert (
        min(
            item["initial_min_nonattached_bead_distance_m"] for item in records.values()
        )
        > 0
    )


def test_neutral_placement_rejects_unachievable_bead_clearance() -> None:
    body = np.array([[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    flag = np.array([[0.0, -0.25, 0.0], [1.0, -0.25, 0.0]])
    with pytest.raises(
        ValueError, match="No collision-free initial hook-neutral placement"
    ):
        _neutralize_initial_hooks_without_bead_overlap(
            body, [flag], np.array([0]), hook_length_um=0.25, bead_diameter_um=10.0
        )
