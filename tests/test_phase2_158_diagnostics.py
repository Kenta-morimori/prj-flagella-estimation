from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pytest

from sim_swim.analysis.flagella_count_behavior import save_state_archive
from sim_swim.analysis.phase2_158_diagnostics import (
    Phase2158DiagnosticConfig,
    _classify_hypotheses,
    _distance_for_row,
    analyze_phase2_158_diagnostics,
)
from sim_swim.sim.core import SimulationState


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _write_fixture(root: Path) -> Path:
    campaign_path = root / "campaign.yaml"
    campaign_path.write_text(
        "\n".join(
            [
                "kind: generic_multi_run",
                "base_config: conf/sim_swim.yaml",
                "base_overrides:",
                "  time.duration_s: 0.003",
                "  time.dt_star: 1.0e-4",
                "sweep:",
                "  axes:",
                "    attach_seed:",
                "      key: seed.attach_seed",
                "      short_name: as",
                "      values: [0]",
                "      ids: [as000]",
                "    phase_seed:",
                "      key: seed.phase_seed",
                "      short_name: ps",
                "      values: [0, 1]",
                "      ids: [ps000, ps001]",
                "    n_flagella:",
                "      key: flagella.n_flagella",
                "      short_name: nf",
                "      values: [3]",
                "      ids: [nf03]",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    input_dir = root / "input"
    summary_rows: list[dict[str, object]] = []
    for phase_seed in (0, 1):
        condition_id = f"as000__ps{phase_seed:03d}__nf03"
        condition_dir = input_dir / condition_id
        states = [
            SimulationState(
                t=index * 0.001,
                position_um=(0.0, 0.0, 0.0),
                quaternion=(0.0, 0.0, 0.0, 1.0),
                velocity_um_s=(0.1, 0.0, 0.0),
                omega_rad_s=(0.0, 0.0, 0.1),
                bead_positions_um=np.zeros((80, 3), dtype=float),
            )
            for index in range(4)
        ]
        save_state_archive(condition_dir / "state_archive.npz", states)
        rows: list[dict[str, object]] = []
        for index, state in enumerate(states):
            fail = phase_seed == 1 and index >= 2
            rows.append(
                {
                    "t_s": state.t,
                    "finite_pass": True,
                    "shape_pass_nonbody_strict": not fail,
                    "first_fail_category_nonbody_strict": "flag" if fail else "none",
                    "body_speed_um_s": 1.0 + index,
                    "body_angular_velocity_rms_rad_s": 2.0 + index,
                    "body_roll_rate_hz": 0.1,
                    "axis_center_body_relative_phase_deg": 10.0 * index,
                    "flag_bond_rel_err_max": 0.2 + index * 0.3,
                    "flag_bond_rel_err_max_flag_id": 1,
                    "flag_bond_rel_err_max_local_bead_i": 1,
                    "flag_bond_rel_err_max_local_bead_j": 2,
                    "flag_bond_rel_err_local_0_1_per_flag": "0:0.1|1:0.2|2:0.3",
                    "flag_bond_rel_err_local_1_2_per_flag": f"0:0.1|1:{0.2 + index * 0.4}|2:0.3",
                    "flag_bond_rel_err_local_2_3_per_flag": "0:0.1|1:0.2|2:0.3",
                    "flag_flag_bead_pair_dist_min_um": 0.5,
                    "flag_flag_close_pair_count": index,
                    "flag_flag_repulsion_force_max_N": 1.0e-14 * index,
                    "motor_torque_Nm": 2.0e-20,
                    "F_motor_mean_body": 1.0e-12,
                    "F_motor_mean_flag": 2.0e-12,
                    "F_spring_mean_flag": 3.0e-12,
                    "F_hook_mean_flag": 4.0e-12,
                    "F_repulsion_mean_flag": 5.0e-12,
                }
            )
        _write_csv(condition_dir / "step_summary.csv", rows)
        summary_rows.append(
            {
                "condition_id": condition_id,
                "body_roll_net_abs_revolutions": 0.1,
                "axis_center_body_relative_net_abs_revolutions_mean": 2.0,
                "axis_center_to_body_roll_ratio_mean": 20.0,
            }
        )
    _write_csv(input_dir / "summary.csv", summary_rows)
    return campaign_path


def test_phase2_158_distance_for_row_prefers_post_step_time() -> None:
    distances = {
        0.0: (10.0, 9.8),
        0.1: (1.0, 0.8),
    }

    assert _distance_for_row({"t_s": "0.0", "dt_s": "0.1"}, distances) == (
        1.0,
        0.8,
    )


def test_phase2_158_contact_hypothesis_reports_measured_counts() -> None:
    run_rows = [
        {
            "n_flagella": "3",
            "strict_pass": "False",
            "attach_seed": "1",
            "first_fail_local_bond": "1-2",
            "first_fail_flag_id": "0",
        },
        {
            "n_flagella": "3",
            "strict_pass": "False",
            "attach_seed": "2",
            "first_fail_local_bond": "2-3",
            "first_fail_flag_id": "1",
        },
        {
            "n_flagella": "3",
            "strict_pass": "True",
            "attach_seed": "0",
            "first_fail_local_bond": "",
            "first_fail_flag_id": "",
        },
    ]
    lead_lag_rows = [
        {
            "window_s": "0.05",
            "interpretation": "bond_growth_without_contact_precursor",
        },
        {
            "window_s": "0.05",
            "interpretation": "bond_growth_with_contact_correlation",
        },
    ]

    assessments = _classify_hypotheses(run_rows, [], lead_lag_rows)
    contact = next(
        item
        for item in assessments
        if item["hypothesis"]
        == "flag-flag or flag-body contact contributes to load variation"
    )

    assert "1/2 failures" in contact["evidence"]
    assert "without a flag-flag or flag-body contact precursor" in contact["evidence"]
    assert contact["bond_growth_without_contact_precursor_count"] == 1
    assert contact["bond_growth_with_contact_correlation_count"] == 1


@pytest.mark.light
def test_phase2_158_diagnostics_writes_summary_events_and_manifest(
    tmp_path: Path,
) -> None:
    campaign_path = _write_fixture(tmp_path)
    output_dir = tmp_path / "out"

    result = analyze_phase2_158_diagnostics(
        Phase2158DiagnosticConfig(
            campaign_config=campaign_path,
            input_dir=tmp_path / "input",
            output_dir=output_dir,
            fail_condition_ids=("as000__ps001__nf03",),
        )
    )

    assert result == output_dir
    summary_rows = list(
        csv.DictReader((output_dir / "run_diagnostic_summary.csv").open())
    )
    assert len(summary_rows) == 2
    failed = [row for row in summary_rows if row["strict_pass"] == "False"]
    assert len(failed) == 1
    assert failed[0]["first_fail_t_s"] == "0.002"
    assert failed[0]["first_fail_local_bond"] == "1-2"
    assert failed[0]["max_proximal_local_bond"] == "1-2"

    event_rows = list(csv.DictReader((output_dir / "failure_event_table.csv").open()))
    assert event_rows
    assert {row["condition_id"] for row in event_rows} == {"as000__ps001__nf03"}
    lead_lag_rows = list(
        csv.DictReader((output_dir / "failure_lead_lag_summary.csv").open())
    )
    assert {row["window_s"] for row in lead_lag_rows} == {"0.25", "0.1", "0.05"}
    assert {
        row["interpretation"] for row in lead_lag_rows if row["window_s"] == "0.05"
    } == {"bond_growth_with_contact_correlation"}
    seed_rows = list(csv.DictReader((output_dir / "seed_failure_table.csv").open()))
    assert len(seed_rows) == 2
    assert {row["condition_id"] for row in seed_rows} == {
        "as000__ps000__nf03",
        "as000__ps001__nf03",
    }
    attach_rows = list(
        csv.DictReader((output_dir / "attach_geometry_table.csv").open())
    )
    assert len(attach_rows) == 1
    assert attach_rows[0]["attach_body_indices"]
    assert attach_rows[0]["attach_pair_dist_min_um"]
    assert (output_dir / "plots" / "as000__ps001__nf03_first_fail_sync.png").is_file()

    manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
    assert manifest["long_simulation_rerun"] is False
    assert manifest["run_counts"] == {"total": 2, "strict_pass": 1, "fail": 1}
    assert manifest["outputs"]["plots"]
    assert manifest["outputs"]["failure_lead_lag_summary"].endswith(
        "failure_lead_lag_summary.csv"
    )
    assert manifest["outputs"]["attach_geometry_table"].endswith(
        "attach_geometry_table.csv"
    )
