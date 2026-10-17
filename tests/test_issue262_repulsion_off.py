"""Issue #262's OFF campaigns must remain paired with the executed #255 model."""

import csv
import json
from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest

from sim_swim.analysis.body_axis_diagnostic import axis_rows, compare_axis_outputs
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
from sim_swim.dynamics.engine import DynamicsEngine
from sim_swim.model.builder import ModelBuilder
from sim_swim.sim.params import SimulationConfig

ROOT = Path(__file__).resolve().parents[1]
EXPECTED_IDS = {
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


@pytest.mark.parametrize("stage", ["1tau", "2s"])
def test_off_campaign_differs_from_issue255_only_by_repulsion(
    stage: str, tmp_path: Path
) -> None:
    old = normalize_campaign_config(
        load_yaml(
            ROOT
            / f"conf/phase2_multi_run/2010_hex_project_local_reaction_{stage}_issue255.yaml"
        )
    )
    new = normalize_campaign_config(
        load_yaml(
            ROOT
            / f"conf/phase2_multi_run/2010_hex_project_repulsion_off_{stage}_issue262.yaml"
        )
    )
    assert new["base_config"] == old["base_config"]
    expected_base = deepcopy(old["base_overrides"])
    expected_base["potentials"]["spring_spring_repulsion"][
        "flagella_flagella_enabled"
    ] = False
    assert new["base_overrides"] == expected_base
    before = {row["condition_id"]: row for row in build_campaign_conditions(old)}
    after = {row["condition_id"]: row for row in build_campaign_conditions(new)}
    assert set(before) == set(after) == EXPECTED_IDS
    assert set(geometry_preflight(new, list(after.values()))) == EXPECTED_IDS
    for cid in EXPECTED_IDS:
        expected_overrides = deepcopy(before[cid]["config_overrides"])
        expected_overrides["potentials"]["spring_spring_repulsion"][
            "flagella_flagella_enabled"
        ] = False
        assert after[cid]["config_overrides"] == expected_overrides
        assert after[cid]["axis_values"] == before[cid]["axis_values"]

    job = load_parallel_job(
        ROOT / f"conf/phase2_parallel/issue262_repulsion_off/hex_{stage}_job.yaml"
    )
    execution = resolve_execution(job, None)
    plan = build_plan(job, execution, tmp_path / "plan")
    assert set(job.condition_ids) == EXPECTED_IDS
    assert job.preflight == "geometry_all_conditions"
    assert execution.worker_policy == "cs10_qualified"
    assert execution.max_workers == 8
    assert len({row["output_dir"] for row in plan["configs"]}) == 13

    base = load_yaml(ROOT / new["base_config"])
    for cid, row in after.items():
        cfg = SimulationConfig.from_dict(base).with_overrides(row["config_overrides"])
        assert cfg.potentials.torsion.force_method == "finite_difference"
        assert cfg.motor.body_reaction_support == "attach_one_ring"
        assert not cfg.potentials.spring_spring_repulsion.body_flagella_enabled
        assert not cfg.potentials.spring_spring_repulsion.flagella_flagella_enabled
        counts = DynamicsEngine(
            ModelBuilder(cfg).build(), cfg
        ).segment_repulsion_pair_counts
        assert counts["body_flagella_active"] == 0
        assert counts["flagella_flagella_active"] == 0
        assert counts["active_total"] == 0


def test_axis_windows_keep_short_direction_excursion_visible() -> None:
    times = np.arange(0.0, 0.301, 0.001)
    angles = np.zeros_like(times)
    angles[(times >= 0.1) & (times < 0.2)] = np.deg2rad(30)
    axes = np.column_stack((np.cos(angles), np.sin(angles), np.zeros_like(times)))
    series, windows, summary = axis_rows(
        condition_id="synthetic", times=times, axes=axes, wobble_final_deg=4.0
    )
    assert len(series) == len(times)
    assert len(windows) == 3
    assert max(row["max_angle_from_initial_deg"] for row in windows) == pytest.approx(
        30
    )
    assert windows[0]["max_angle_from_initial_deg"] == 0
    assert summary["body_axis_wobble_rms_deg_final"] == 4.0


def test_axis_windows_ignore_integrator_endpoint_roundoff() -> None:
    times = np.r_[np.arange(0, 2, 0.001), 2.00000000002]
    axes = np.tile([1.0, 0.0, 0.0], (len(times), 1))
    _, windows, _ = axis_rows(
        condition_id="synthetic", times=times, axes=axes, wobble_final_deg=0.0
    )
    assert len(windows) == 20
    assert windows[-1]["window_end_s"] == 2.0


def test_pair_windows_requires_matching_config_and_condition_coverage(
    tmp_path: Path,
) -> None:
    old_config = Path(
        "conf/phase2_multi_run/2010_hex_project_local_reaction_2s_issue255.yaml"
    )
    off_config = Path(
        "conf/phase2_multi_run/2010_hex_project_repulsion_off_2s_issue262.yaml"
    )
    old_dir, off_dir = tmp_path / "old", tmp_path / "off"
    for directory, config, value in (
        (old_dir, old_config, 1.0),
        (off_dir, off_config, 2.0),
    ):
        directory.mkdir()
        (directory / "manifest.json").write_text(
            json.dumps(
                {
                    "condition_count": 13,
                    "source_config": str(config),
                }
            )
        )
        with (directory / "axis_windows.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "condition_id",
                    "window_start_s",
                    "window_end_s",
                    "direction_rms_deg",
                    "direction_max_deg",
                    "endpoint_drift_deg",
                ],
            )
            writer.writeheader()
            for cid in sorted(EXPECTED_IDS):
                writer.writerow(
                    {
                        "condition_id": cid,
                        "window_start_s": 0.0,
                        "window_end_s": 0.1,
                        "direction_rms_deg": value,
                        "direction_max_deg": value,
                        "endpoint_drift_deg": value,
                    }
                )
    result = compare_axis_outputs(
        reference_config=old_config,
        off_config=off_config,
        reference_axis_dir=old_dir,
        off_axis_dir=off_dir,
        output_dir=tmp_path / "pair",
    )
    with result.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 13
    assert all(float(row["off_minus_on_direction_rms_deg"]) == 1.0 for row in rows)
    with (off_dir / "axis_windows.csv").open("a") as handle:
        handle.write("unexpected,0,0.1,2,2,2\n")
    with pytest.raises(ValueError, match="time windows do not match"):
        compare_axis_outputs(
            reference_config=old_config,
            off_config=off_config,
            reference_axis_dir=old_dir,
            off_axis_dir=off_dir,
            output_dir=tmp_path / "pair",
        )
