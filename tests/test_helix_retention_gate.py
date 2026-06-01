from __future__ import annotations

from sim_swim.sim.helix_retention_gate import (
    summarize_single_flagellum_helix_retention,
)


def _base_row(
    step: int,
    phase_deg: float,
    *,
    root_phase_deg: float | None = None,
) -> dict[str, str]:
    root_phase = phase_deg if root_phase_deg is None else root_phase_deg
    return {
        "step": str(step),
        "t_s": f"{step * 0.001:.6f}",
        "finite_pass": "True",
        "shape_pass_nonbody": "True",
        "first_fail_category_nonbody": "none",
        "flag_phase_deg": f"{root_phase:.6f}",
        "flag_phase_rate_hz": "0.0",
        "flag_helix_spin_phase_deg": f"{phase_deg:.6f}",
        "flag_helix_spin_rate_hz": "10.0",
        "flag_helix_spin_fit_r2": "0.99",
        "hook_len_rel_err_max": "0.01",
        "local_attach_first_rel_err": "0.01",
        "flag_bond_rel_err_max": "0.01",
        "flag_bend_err_max_deg": "1.0",
        "flag_torsion_err_max_deg": "1.0",
        "motor_degenerate_axis_count": "0",
        "motor_bond_length_clipped_count": "0",
    }


def test_helix_retention_gate_requires_net_spin_not_only_instant_rate() -> None:
    phases = [0.0, 20.0, -20.0, 20.0, -20.0, 0.0] * 20
    rows = [_base_row(step, phase) for step, phase in enumerate(phases)]

    summary = summarize_single_flagellum_helix_retention(
        rows,
        skip_initial_steps=0,
        min_steps=50,
        min_net_abs_spin_revolutions=0.5,
    )

    assert summary["helix_retention_pass"] is False
    assert summary["first_fail_category"] == "motor_no_rotation"
    assert summary["median_abs_flag_helix_spin_rate_hz"] == 10.0
    assert summary["net_abs_flag_root_revolutions"] == 0.0
    assert summary["net_abs_flag_helix_spin_revolutions"] == 0.0
    assert summary["flag_helix_spin_direction_consistency"] == 0.0


def test_helix_retention_gate_accepts_monotonic_net_spin() -> None:
    rows = [_base_row(step, step * 10.0) for step in range(120)]

    summary = summarize_single_flagellum_helix_retention(
        rows,
        skip_initial_steps=0,
        min_steps=50,
        min_net_abs_spin_revolutions=1.0,
    )

    assert summary["helix_retention_pass"] is True
    assert summary["first_fail_category"] == "none"
    assert summary["net_abs_flag_helix_spin_revolutions"] > 3.0
    assert summary["flag_helix_spin_direction_consistency"] == 1.0
    assert summary["helix_to_root_net_rotation_ratio"] == 1.0


def test_helix_retention_summary_reports_root_to_helix_transfer_ratio() -> None:
    rows = [
        _base_row(step, step * 1.0, root_phase_deg=step * 10.0) for step in range(120)
    ]

    summary = summarize_single_flagellum_helix_retention(
        rows,
        skip_initial_steps=0,
        min_steps=50,
        min_net_abs_spin_revolutions=0.1,
    )

    assert summary["helix_retention_pass"] is True
    assert summary["net_abs_flag_root_revolutions"] > 3.0
    assert summary["net_abs_flag_helix_spin_revolutions"] > 0.3
    assert 0.09 < summary["helix_to_root_net_rotation_ratio"] < 0.11
