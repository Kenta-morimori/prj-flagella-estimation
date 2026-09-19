from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from sim_swim.analysis.model_development_evaluation import (
    build_evaluation,
    collect_rows,
)
from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    load_yaml,
    normalize_campaign_config,
)


ROOT = Path(__file__).resolve().parents[1]
CONFIG = ROOT / "conf/phase2_multi_run/2010_hex_project_torque_dt_1tau_issue244.yaml"


def _summary(*, hook_angle_only: bool = False, hook_length_fail: bool = False) -> dict:
    metrics = {
        "hook_angle_err_max_deg": {"max": 31.0 if hook_angle_only else 1.0},
        "local_attach_first_rel_err": {"max": 0.01},
        "hook_len_rel_err_max": {"max": 1.1 if hook_length_fail else 0.01},
        "flag_bond_rel_err_max": {"max": 0.03},
        "flag_bend_err_max_deg": {"max": 1.0},
        "flag_torsion_err_max_deg": {"max": 1.0},
        "body_spring_max_stretch_ratio": {"max": 0.04},
        "motor_force_balance_residual_ratio": {"max": 1.0e-12},
        "motor_torque_balance_residual_ratio": {"max": 0.2},
    }
    return {
        "gates": {
            "finite": {"any_fail": False},
            "shape_body": {"any_fail": False},
            "shape_nonbody": {
                "any_fail": hook_angle_only or hook_length_fail,
                "first_failure_category": (
                    "hook" if hook_angle_only or hook_length_fail else None
                ),
            },
        },
        "all_step_metrics": metrics,
    }


def _write_runs(root: Path) -> tuple[Path, Path]:
    config = normalize_campaign_config(load_yaml(CONFIG))
    conditions = build_campaign_conditions(config)
    result: list[Path] = []
    for dt_star, name in ((1.0e-4, "fine"), (1.0e-3, "coarse")):
        run_dir = root / name
        run_dir.mkdir()
        records = []
        csv_rows = []
        for condition in conditions:
            if condition["axis_values"]["dt_star"] != dt_star:
                continue
            condition_dir = run_dir / condition["condition_id"]
            condition_dir.mkdir()
            (condition_dir / "run_summary.json").write_text(
                json.dumps(_summary(hook_angle_only=dt_star == 1.0e-4)),
                encoding="utf-8",
            )
            (condition_dir / "performance.json").write_text(
                json.dumps({"wall_time_s": 10.0, "steps_per_s": 100.0}),
                encoding="utf-8",
            )
            records.append(
                {
                    **condition,
                    "output_dir": str(condition_dir),
                    "source_config_path": "conf/sim_swim_2010_hex.yaml",
                    "time": {"dt_star": dt_star, "dt_internal_s": dt_star * 0.04},
                }
            )
            csv_rows.append(
                {"condition_id": condition["condition_id"], "completed": "True"}
            )
        with (run_dir / "summary.csv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.DictWriter(handle, fieldnames=["condition_id", "completed"])
            writer.writeheader()
            writer.writerows(csv_rows)
        (run_dir / "run_manifest.json").write_text(
            json.dumps(
                {
                    "base_config": "conf/sim_swim_2010_hex.yaml",
                    "model_profile": load_yaml(ROOT / "conf/sim_swim_2010_hex.yaml")[
                        "model_profile"
                    ],
                    "git": {"commit": name},
                    "conditions": records,
                }
            ),
            encoding="utf-8",
        )
        result.append(run_dir)
    return result[0], result[1]


def test_pending_profiles_declare_development_evaluation_contract() -> None:
    for path in (
        ROOT / "conf/sim_swim_2010_hex.yaml",
        ROOT / "conf/sim_swim_2015.yaml",
        ROOT / "conf/sim_swim_2015_paper.yaml",
    ):
        raw = load_yaml(path)
        assert raw["model_profile"]["implementation_status"] == "pending"
        assert raw["development_evaluation"]["version"] == 1
        assert raw["development_evaluation"]["feature_evaluation"] == "separate_issue"


def test_issue244_screen_is_60_conditions() -> None:
    config = normalize_campaign_config(load_yaml(CONFIG))
    conditions = build_campaign_conditions(config)
    assert len(conditions) == 60
    assert conditions[0]["condition_id"] == "nf01__tq1p0e20__dt1e4"
    assert conditions[-1]["condition_id"] == "nf06__tq3p5e20__dt1e3"
    assert config["development_evaluation"]["expected_condition_count"] == 60


def test_common_evaluation_reuses_runs_and_ignores_hook_angle_only(
    tmp_path: Path,
) -> None:
    fine, coarse = _write_runs(tmp_path)
    rows, provenance = collect_rows(config=load_yaml(CONFIG), run_dirs=[fine, coarse])
    assert len(rows) == 60
    assert {row["screen_status"] for row in rows} == {"pass"}
    assert {row["raw_nonbody_any_fail"] for row in rows} == {False, True}
    assert len(provenance) == 2
    outputs = build_evaluation(
        config_path=CONFIG,
        run_dirs=[fine, coarse],
        output_dir=tmp_path / "model_development_evaluation",
    )
    assert outputs["summary_csv"].is_file()
    assert outputs["manifest"].is_file()
    for n_flagella in range(1, 7):
        assert outputs[f"heatmap_nf{n_flagella:02d}"].is_file()
    manifest = json.loads(outputs["manifest"].read_text(encoding="utf-8"))
    assert manifest["status_counts"] == {"pass": 60, "fail": 0}


def test_common_evaluation_rejects_remaining_hook_length_failure(
    tmp_path: Path,
) -> None:
    fine, coarse = _write_runs(tmp_path)
    target = fine / "nf01__tq1p0e20__dt1e4" / "run_summary.json"
    target.write_text(json.dumps(_summary(hook_length_fail=True)), encoding="utf-8")
    rows, _ = collect_rows(config=load_yaml(CONFIG), run_dirs=[fine, coarse])
    selected = next(
        row for row in rows if row["condition_id"] == "nf01__tq1p0e20__dt1e4"
    )
    assert selected["screen_status"] == "fail"


def test_common_evaluation_rejects_missing_cell(tmp_path: Path) -> None:
    fine, coarse = _write_runs(tmp_path)
    manifest_path = coarse / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["conditions"] = manifest["conditions"][1:]
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="Missing expected conditions"):
        collect_rows(config=load_yaml(CONFIG), run_dirs=[fine, coarse])
