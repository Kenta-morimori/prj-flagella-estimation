import json
from pathlib import Path

import numpy as np
import pytest

from sim_swim.analysis.multi_run_campaign import load_yaml
from sim_swim.analysis.torque_dt_stability import build_plan
from sim_swim.analysis.torque_dt_stability_campaign import (
    _assert_condition,
    _validate_campaign_contract,
    summarize_campaign,
)
from sim_swim.analysis.sweeps.generic_multi_run import run_campaign
from sim_swim.sim.core import Simulator
from sim_swim.sim.params import SimulationConfig


CONFIG = Path("conf/phase2_multi_run/2010_project_torque_dt_initial_screen.yaml")


def _summary() -> dict[str, object]:
    return {
        "execution": {"status": "completed"},
        "gates": {
            name: {"any_fail": False}
            for name in ("finite", "shape_nonbody", "shape_body")
        },
        "all_step_metrics": {
            "motor_force_balance_residual_ratio": {"max": 1.0e-12},
            "motor_torque_balance_residual_ratio": {"max": 1.0e-12},
        },
    }


def _campaign_fixture(root: Path) -> None:
    raw = load_yaml(CONFIG)
    base = load_yaml(Path(raw["base_config"]))
    plan = build_plan(CONFIG)
    records = []
    for condition in plan["conditions"]:
        cfg = SimulationConfig.from_dict(base).with_overrides(
            condition["config_overrides"]
        )
        directory = root / condition["condition_id"]
        directory.mkdir(parents=True)
        model = Simulator(cfg).model
        frame = model.positions_m * 1.0e6
        samples = int(condition["comparison_sample_count"])
        np.savez_compressed(
            directory / "state_archive.npz",
            t=np.linspace(0.0, cfg.tau_s, samples),
            position_um=np.zeros((samples, 3)),
            quaternion=np.tile(np.array([0.0, 0.0, 0.0, 1.0]), (samples, 1)),
            velocity_um_s=np.zeros((samples, 3)),
            bead_positions_um=np.tile(frame, (samples, 1, 1)),
        )
        (directory / "run_summary.json").write_text(json.dumps(_summary()))
        records.append(
            {
                "condition_id": condition["condition_id"],
                "torque_Nm_per_flagellum": condition["torque_Nm_per_flagellum"],
                "dt_star": condition["dt_star"],
                "comparison_role": condition["comparison_role"],
                "output_dir": str(directory),
                "config_overrides": condition["config_overrides"],
                "comparison_sample_count": samples,
            }
        )
    (root / "run_manifest.json").write_text(json.dumps({"conditions": records}))


def test_initial_screen_is_diagnostic_only_and_records_similarity(
    tmp_path: Path,
) -> None:
    _campaign_fixture(tmp_path)
    payload = summarize_campaign(tmp_path, config_path=CONFIG)

    assert len(payload["conditions"]) == 8
    assert len(payload["dt_comparisons"]) == 4
    assert {row["status"] for row in payload["dt_comparisons"]} == {"diagnostic_only"}
    assert payload["interpretation"]["formal_dt_adoption_status"] == (
        "awaiting_1e-5_reference"
    )
    assert len(payload["torque_similarity"]) == 3
    assert (tmp_path / "summary.csv").is_file()
    assert (tmp_path / "dt_comparison.csv").is_file()
    assert (tmp_path / "torque_similarity.csv").is_file()
    assert (tmp_path / "qc_summary.json").is_file()


def test_campaign_marks_bad_reference_not_comparable(tmp_path: Path) -> None:
    _campaign_fixture(tmp_path)
    manifest_path = tmp_path / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    reference = next(
        row
        for row in manifest["conditions"]
        if row["comparison_role"] == "screen_comparator"
    )
    summary_path = Path(reference["output_dir"]) / "run_summary.json"
    failed = _summary()
    failed["gates"]["finite"]["any_fail"] = True  # type: ignore[index]
    summary_path.write_text(json.dumps(failed))
    payload = summarize_campaign(tmp_path, config_path=CONFIG)

    affected = [
        row
        for row in payload["dt_comparisons"]
        if row["reference_condition_id"] == reference["condition_id"]
    ]
    assert {row["status"] for row in affected} == {"not_comparable"}


def test_contract_rejects_torque_mismatch_and_generic_dry_run_starts_no_simulation() -> (
    None
):
    raw = load_yaml(CONFIG)
    base = load_yaml(Path(raw["base_config"]))
    condition = build_plan(CONFIG)["conditions"][0]
    cfg = SimulationConfig.from_dict(base).with_overrides(condition["config_overrides"])
    _assert_condition(cfg, condition)
    bad = cfg.with_overrides({"motor": {"reference_torque_Nm": 2.0e-20}})
    with pytest.raises(ValueError, match="torques"):
        _assert_condition(bad, condition)
    _validate_campaign_contract(raw)
    assert run_campaign(["--campaign-config", str(CONFIG), "--dry-run"]) == Path()


def test_initial_screen_rejects_partial_execution_and_cli_overrides() -> None:
    with pytest.raises(ValueError, match="fixed condition contract"):
        run_campaign(["--campaign-config", str(CONFIG), "--sample-limit", "1"])
    with pytest.raises(ValueError, match="fixed condition contract"):
        run_campaign(
            [
                "--campaign-config",
                str(CONFIG),
                "time.integration.dt_star=1e-5",
            ]
        )
