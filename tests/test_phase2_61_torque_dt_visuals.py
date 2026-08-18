from __future__ import annotations

import json
from pathlib import Path

import pytest

from sim_swim.analysis.torque_dt_stability_visuals import build_visuals, feature_rows


pytestmark = pytest.mark.light


def _summary() -> dict[str, object]:
    return {
        "execution": {"status": "completed"},
        "gates": {
            "finite": {"any_fail": False},
            "shape_nonbody": {"any_fail": False},
            "shape_body": {"any_fail": False},
        },
        "all_step_metrics": {
            "body_spring_max_stretch_ratio": {"max": 0.01},
            "flag_bond_rel_err_max": {"max": 0.01},
            "hook_len_rel_err_max": {"max": 0.02},
            "flag_helix_axis_mean_deviation_deg_max": {"max": 3.0},
            "flag_helix_axis_rearward_projection_min": {"min": 0.99},
            "bundle_participation_ratio": {"final": 2.0 / 3.0},
            "flag_helix_bundle_radius_max_um": {"max": 1.5},
        },
    }


def _fixture(run_dir: Path) -> None:
    conditions = []
    for torque in (1.0e-21, 1.2e-18):
        for dt_star in (1.0e-3, 1.0e-4):
            condition_id = f"T{torque:.0e}_dt{dt_star:.0e}"
            condition_dir = run_dir / condition_id
            condition_dir.mkdir(parents=True)
            (condition_dir / "run_summary.json").write_text(
                json.dumps(_summary()), encoding="utf-8"
            )
            conditions.append(
                {
                    "condition_id": condition_id,
                    "output_dir": str(condition_dir),
                    "torque_Nm_per_flagellum": torque,
                    "dt_star": dt_star,
                    "time": {"tau_s": 1.0e-3, "duration_s": 1.0e-3},
                }
            )
    (run_dir / "run_manifest.json").write_text(
        json.dumps({"conditions": conditions}), encoding="utf-8"
    )
    (run_dir / "qc_summary.json").write_text(
        json.dumps(
            {
                "conditions": [
                    {"condition_id": condition["condition_id"], "status": "pass"}
                    for condition in conditions
                ]
            }
        ),
        encoding="utf-8",
    )


def test_visuals_extract_features_and_write_heatmap(tmp_path: Path) -> None:
    _fixture(tmp_path)

    rows = feature_rows(tmp_path)
    assert len(rows) == 4
    assert all(row["required_qc_pass"] == 1 for row in rows)
    assert rows[0]["max_axis_mean_deviation_deg"] == pytest.approx(3.0)

    outputs = build_visuals(tmp_path, tmp_path / "analysis")
    assert outputs["csv"].is_file()
    assert outputs["heatmap"].is_file()
    assert outputs["manifest"].is_file()


def test_visuals_split_issue199_heatmaps_by_tau_policy(tmp_path: Path) -> None:
    _fixture(tmp_path)
    manifest_path = tmp_path / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for index, condition in enumerate(manifest["conditions"]):
        condition["tau_policy"] = (
            "tau_fixed_control" if index % 2 == 0 else "torque_linked_tau"
        )
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    outputs = build_visuals(tmp_path, tmp_path / "analysis")
    payload = json.loads(outputs["manifest"].read_text(encoding="utf-8"))

    assert len(payload["outputs"]["feature_heatmaps_png"]) == 2
    assert all(
        Path(path).is_file() for path in payload["outputs"]["feature_heatmaps_png"]
    )
