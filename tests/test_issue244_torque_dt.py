from __future__ import annotations

import json
from pathlib import Path

from sim_swim.analysis.issue244_torque_dt import build_analysis, feature_rows
from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    load_yaml,
    normalize_campaign_config,
)


ROOT = Path(__file__).resolve().parents[1]


def _summary(*, warning: bool = False, persistent_hook: bool = False) -> dict:
    nonbody = {
        "final_pass": True,
        "any_fail": warning or persistent_hook,
        "first_observed_fail_t_s": 4.0e-6,
        "observed_fail_sample_count": 2 if persistent_hook else (1 if warning else 0),
        "first_failure_category": "hook" if warning or persistent_hook else None,
        "first_failure_snapshot": {
            "step": 0,
            "metrics": {
                "hook_angle_err_max_deg": 31.0,
                "local_attach_first_rel_err": 0.01,
                "hook_len_rel_err_max": 0.01,
            },
        }
        if warning or persistent_hook
        else None,
    }
    metrics = {
        "hook_angle_err_max_deg": {"max": 31.0},
        "flag_bond_rel_err_max": {"max": 0.03},
        "body_spring_max_stretch_ratio": {"max": 0.04},
        "motor_force_balance_residual_ratio": {"max": 1.0e-12},
        "motor_torque_balance_residual_ratio": {"max": 0.2},
    }
    return {
        "gates": {
            "finite": {"any_fail": False},
            "shape_nonbody": nonbody,
            "shape_body": {"any_fail": False},
        },
        "all_step_metrics": metrics,
    }


def _write_run(root: Path, *, dt_star: float, warning: bool = False) -> None:
    records = []
    for n_flagella in (1, 4):
        for torque in (1.0e-20, 2.0e-20, 2.5e-20, 3.0e-20, 3.5e-20):
            condition = root / f"nf{n_flagella:02d}_{torque:.1e}"
            condition.mkdir(parents=True)
            (condition / "run_summary.json").write_text(
                json.dumps(_summary(warning=warning)), encoding="utf-8"
            )
            (condition / "performance.json").write_text(
                json.dumps({"wall_time_s": 10.0, "steps_per_s": 100.0}),
                encoding="utf-8",
            )
            records.append(
                {
                    "condition_id": condition.name,
                    "output_dir": str(condition),
                    "axis_values": {"n_flagella": n_flagella, "motor_torque": torque},
                    "time": {"dt_star": dt_star, "dt_internal_s": dt_star * 0.04},
                }
            )
    (root / "run_manifest.json").write_text(
        json.dumps({"conditions": records}), encoding="utf-8"
    )


def test_issue244_torque_dt_campaign_and_followup_grids_are_complete() -> None:
    screen = normalize_campaign_config(
        load_yaml(
            ROOT / "conf/phase2_multi_run/2010_hex_project_torque_dt_1tau_issue244.yaml"
        )
    )
    seed_grid = normalize_campaign_config(
        load_yaml(
            ROOT / "conf/phase2_multi_run/2010_hex_project_seed_grid_1tau_issue244.yaml"
        )
    )
    convergence = normalize_campaign_config(
        load_yaml(
            ROOT
            / "conf/phase2_multi_run/2010_hex_project_dt_convergence_1tau_issue244.yaml"
        )
    )
    conditions = build_campaign_conditions(screen)
    assert len(conditions) == 20
    assert conditions[0]["condition_id"] == "nf01__tq1p0e20__dt1e4"
    assert conditions[-1]["condition_id"] == "nf04__tq3p5e20__dt1e3"
    assert len(build_campaign_conditions(seed_grid)) == 54
    assert len(build_campaign_conditions(convergence)) == 4
    assert {
        condition["config_overrides"]["time"]["integration"]["dt_star"]
        for condition in build_campaign_conditions(convergence)
    } == {1.0e-4, 5.0e-5}


def test_issue244_analysis_writes_two_heatmaps_and_initial_hook_warning(
    tmp_path: Path,
) -> None:
    baseline = tmp_path / "baseline"
    coarse = tmp_path / "coarse"
    _write_run(baseline, dt_star=1.0e-4, warning=True)
    _write_run(coarse, dt_star=1.0e-3, warning=False)
    rows = feature_rows(baseline_run_dir=baseline, coarse_run_dir=coarse)
    assert len(rows) == 20
    assert {row["screen_status"] for row in rows} == {"warning", "pass"}
    outputs = build_analysis(
        baseline_run_dir=baseline,
        coarse_run_dir=coarse,
        output_dir=tmp_path / "analysis",
    )
    assert outputs["summary_csv"].is_file()
    assert outputs["heatmap_nf01"].is_file()
    assert outputs["heatmap_nf04"].is_file()
    assert outputs["manifest"].is_file()


def test_issue244_analysis_rejects_persistent_hook_failure(tmp_path: Path) -> None:
    baseline = tmp_path / "baseline"
    coarse = tmp_path / "coarse"
    _write_run(baseline, dt_star=1.0e-4, warning=True)
    _write_run(coarse, dt_star=1.0e-3, warning=False)
    target = coarse / "nf01_1.0e-20" / "run_summary.json"
    target.write_text(json.dumps(_summary(persistent_hook=True)), encoding="utf-8")
    rows = feature_rows(baseline_run_dir=baseline, coarse_run_dir=coarse)
    selected = next(
        row
        for row in rows
        if row["n_flagella"] == 1
        and row["torque_Nm"] == 1.0e-20
        and row["dt_star"] == 1.0e-3
    )
    assert selected["screen_status"] == "fail"
