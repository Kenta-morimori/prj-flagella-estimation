from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import pytest

from sim_swim.analysis.multi_run_campaign import (
    build_campaign_conditions,
    load_yaml,
    normalize_campaign_config,
)
from sim_swim.analysis.spring_formulation_comparison import (
    decide_default,
    force_extension_rows,
    write_force_extension_artifacts,
)
from sim_swim.sim.params import SimulationConfig

ROOT = Path(__file__).resolve().parents[1]


def _write_summary(
    path: Path,
    *,
    fene_pass: bool,
    legacy_pass: bool = True,
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "axis_formulation_value",
                "duration_s",
                "final_t_s",
                "final_shape_pass_nonbody",
                "body_shape_pass",
                "first_fail_t_s",
                "max_flag_bond_rel_err",
                "max_hook_len_rel_err",
            ],
        )
        writer.writeheader()
        for formulation, passed in (
            ("legacy", legacy_pass),
            ("fene_fraenkel", fene_pass),
        ):
            writer.writerow(
                {
                    "axis_formulation_value": formulation,
                    "duration_s": "1.0",
                    "final_t_s": "1.0",
                    "final_shape_pass_nonbody": str(passed),
                    "body_shape_pass": str(passed),
                    "first_fail_t_s": "" if passed else "0.25",
                    "max_flag_bond_rel_err": "0.01",
                    "max_hook_len_rel_err": "0.02",
                }
            )


@pytest.mark.light
@pytest.mark.parametrize(
    ("profile_name", "duration_s", "torque_nm"),
    [
        ("spring_formulation_motor_off.yaml", 0.1, 0.0),
        ("spring_formulation_motor_on.yaml", 1.0, 2.5e-20),
    ],
)
def test_spring_formulation_profiles_expand_two_conditions(
    profile_name: str,
    duration_s: float,
    torque_nm: float,
) -> None:
    profile = normalize_campaign_config(
        load_yaml(ROOT / "conf/phase2_multi_run" / profile_name)
    )
    conditions = build_campaign_conditions(profile)

    assert [condition["condition_id"] for condition in conditions] == [
        "legacy",
        "fene_fraenkel",
    ]
    assert profile["base_overrides"]["time"]["duration_s"] == duration_s
    assert profile["base_overrides"]["time"]["dt_star"] == 1.0e-4
    assert profile["base_overrides"]["motor"]["torque_Nm"] == torque_nm
    assert profile["base_overrides"]["motor"]["enable_switching"] is False
    assert (
        conditions[1]["config_overrides"]["potentials"]["spring"]["formulation"]
        == "fene_fraenkel"
    )
    base_config = load_yaml(ROOT / str(profile["base_config"]))
    effective_configs = [
        SimulationConfig.from_dict(base_config).with_overrides(
            condition["config_overrides"]
        )
        for condition in conditions
    ]
    assert [cfg.potentials.spring.formulation for cfg in effective_configs] == [
        "legacy",
        "fene_fraenkel",
    ]


@pytest.mark.light
def test_force_extension_uses_expected_equilibrium_and_relative_fene_curve() -> None:
    rows = force_extension_rows(samples=199)

    equilibrium = [
        row for row in rows if abs(float(row["relative_extension_q"])) < 1.0e-12
    ]
    assert len(equilibrium) == 4
    assert all(abs(float(row["force_N"])) < 1.0e-30 for row in equilibrium)

    fene_at_q = [
        row
        for row in rows
        if row["formulation"] == "fene_fraenkel"
        and abs(float(row["relative_extension_q"]) - 0.05) < 1.0e-12
    ]
    assert len(fene_at_q) == 2
    assert float(fene_at_q[0]["force_N"]) == pytest.approx(
        float(fene_at_q[1]["force_N"])
    )


@pytest.mark.light
def test_default_decision_prefers_fene_only_when_both_runs_pass(
    tmp_path: Path,
) -> None:
    motor_off = tmp_path / "motor_off.csv"
    motor_on = tmp_path / "motor_on.csv"
    _write_summary(motor_off, fene_pass=True)
    _write_summary(motor_on, fene_pass=True)

    assert decide_default(motor_off, motor_on)["selected_default"] == ("fene_fraenkel")

    _write_summary(motor_on, fene_pass=False)
    decision = decide_default(motor_off, motor_on)
    assert decision["selected_default"] == "legacy"
    assert decision["evaluations"]["fene_fraenkel"]["overall_pass"] is False


@pytest.mark.light
def test_missing_formulation_row_retains_legacy(tmp_path: Path) -> None:
    motor_off = tmp_path / "motor_off.csv"
    motor_on = tmp_path / "motor_on.csv"
    _write_summary(motor_off, fene_pass=True)
    _write_summary(motor_on, fene_pass=True)
    rows = list(csv.DictReader(motor_on.open(encoding="utf-8", newline="")))
    with motor_on.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerow(rows[0])

    decision = decide_default(motor_off, motor_on)
    assert decision["selected_default"] == "legacy"
    assert decision["evaluations"]["fene_fraenkel"]["motor_on"]["reasons"] == [
        "summary row is missing"
    ]


@pytest.mark.light
def test_incomplete_run_retains_legacy(tmp_path: Path) -> None:
    motor_off = tmp_path / "motor_off.csv"
    motor_on = tmp_path / "motor_on.csv"
    _write_summary(motor_off, fene_pass=True)
    _write_summary(motor_on, fene_pass=True)
    rows = list(csv.DictReader(motor_on.open(encoding="utf-8", newline="")))
    rows[1]["final_t_s"] = "0.5"
    with motor_on.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    decision = decide_default(motor_off, motor_on)
    assert decision["selected_default"] == "legacy"
    assert (
        "final_t_s did not reach duration_s"
        in decision["evaluations"]["fene_fraenkel"]["motor_on"]["reasons"]
    )


@pytest.mark.light
def test_non_finite_summary_metric_retains_legacy(tmp_path: Path) -> None:
    motor_off = tmp_path / "motor_off.csv"
    motor_on = tmp_path / "motor_on.csv"
    _write_summary(motor_off, fene_pass=True)
    _write_summary(motor_on, fene_pass=True)
    rows = list(csv.DictReader(motor_on.open(encoding="utf-8", newline="")))
    rows[1]["max_hook_len_rel_err"] = "nan"
    with motor_on.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    decision = decide_default(motor_off, motor_on)

    assert decision["selected_default"] == "legacy"
    assert (
        "max_hook_len_rel_err is non-finite"
        in decision["evaluations"]["fene_fraenkel"]["motor_on"]["reasons"]
    )


@pytest.mark.light
def test_comparison_cli_writes_run_contract(tmp_path: Path, monkeypatch) -> None:
    motor_off = tmp_path / "motor_off.csv"
    motor_on = tmp_path / "motor_on.csv"
    _write_summary(motor_off, fene_pass=True)
    _write_summary(motor_on, fene_pass=True)

    script_path = ROOT / "scripts/01_simulate_swimming/analyze_spring_formulations.py"
    spec = importlib.util.spec_from_file_location(
        "analyze_spring_formulations", script_path
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr("sim_swim.core.run_context._require_clean_git", lambda: None)

    output_root = module.main(
        [
            f"motor_off_summary={motor_off}",
            f"motor_on_summary={motor_on}",
            f"output_base_dir={tmp_path / 'output'}",
            "samples=9",
        ]
    )

    assert (output_root / "run.log").is_file()
    assert (output_root / "force_extension.csv").is_file()
    assert (output_root / "force_extension.png").is_file()
    assert (output_root / "default_decision.json").is_file()
    assert (output_root / "default_decision.md").is_file()
    manifest = json.loads((output_root / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["kind"] == "spring_formulation_comparison"
    assert manifest["run_time"]["timezone"] == "Asia/Tokyo"
    assert manifest["decision"]["selected_default"] == "fene_fraenkel"


@pytest.mark.light
def test_force_extension_writer_creates_csv_and_png(tmp_path: Path) -> None:
    csv_path, png_path = write_force_extension_artifacts(tmp_path, samples=9)
    assert csv_path.is_file()
    assert png_path.is_file()
