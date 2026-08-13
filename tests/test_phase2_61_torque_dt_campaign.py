import json
import logging
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from sim_swim.analysis.multi_run_campaign import load_yaml
from sim_swim.analysis.torque_dt_stability import build_plan
from sim_swim.analysis.torque_dt_stability_campaign import (
    _assert_condition,
    _validate_campaign_contract,
    summarize_campaign,
)
from sim_swim.analysis import torque_dt_stability_campaign as campaign
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


def test_campaign_motor_reaction_is_balanced_in_a_representative_step(
    tmp_path: Path,
) -> None:
    raw = load_yaml(CONFIG)
    base = load_yaml(Path(raw["base_config"]))
    condition = next(
        row
        for row in build_plan(CONFIG)["conditions"]
        if row["condition_id"] == "T2p5e-20_dt1e-4"
    )
    cfg = SimulationConfig.from_dict(base).with_overrides(condition["config_overrides"])

    Simulator(cfg).run(cfg.dt_s, step_summary_dir=tmp_path)

    summary = json.loads((tmp_path / "run_summary.json").read_text())
    metrics = summary["all_step_metrics"]
    assert metrics["motor_force_balance_residual_ratio"]["max"] <= 1.0e-8
    assert metrics["motor_torque_balance_residual_ratio"]["max"] <= 1.0e-8


def test_campaign_persists_condition_boundary_progress(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The live manifest identifies the active condition before it is run."""

    logger = logging.getLogger("test_phase2_61_campaign_progress")
    logger.handlers.clear()
    logger.addHandler(logging.NullHandler())
    active_condition_ids: list[str] = []

    def fake_init_run(*args: object, **kwargs: object) -> SimpleNamespace:
        (tmp_path / "manifest.json").write_text("{}", encoding="utf-8")
        return SimpleNamespace(out=SimpleNamespace(root=tmp_path), logger=logger)

    def fake_run_condition(
        root: Path,
        base: dict[str, object],
        condition: dict[str, object],
        qc: dict[str, object],
    ) -> dict[str, object]:
        del base
        manifest = json.loads((root / "run_manifest.json").read_text(encoding="utf-8"))
        active = [
            row["condition_id"]
            for row in manifest["conditions"]
            if row.get("execution", {}).get("status") == "running"
        ]
        assert active == [condition["condition_id"]]
        active_condition_ids.extend(active)
        return {
            "condition_id": condition["condition_id"],
            "torque_Nm_per_flagellum": condition["torque_Nm_per_flagellum"],
            "dt_star": condition["dt_star"],
            "comparison_role": condition["comparison_role"],
            "output_dir": str(root / str(condition["condition_id"])),
            "config_overrides": condition["config_overrides"],
            "qc": qc,
        }

    monkeypatch.setattr(campaign, "init_run", fake_init_run)
    monkeypatch.setattr(campaign, "_run_condition", fake_run_condition)
    monkeypatch.setattr(campaign, "summarize_campaign", lambda *args, **kwargs: {})

    assert campaign.run_torque_linked_campaign(CONFIG, output_root=tmp_path) == tmp_path

    manifest = json.loads((tmp_path / "run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["execution"]["status"] == "completed"
    assert len(active_condition_ids) == 8
    for record in manifest["conditions"]:
        execution = record["execution"]
        assert execution["status"] == "completed"
        assert execution["started_at_jst"].endswith("+09:00")
        assert execution["finished_at_jst"].endswith("+09:00")
        assert execution["wall_seconds"] >= 0.0
