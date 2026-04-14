from __future__ import annotations

import pytest

from sim_swim.sim.params import SimulationConfig


def _base_cfg() -> dict:
    return {
        "scale": {"b_um": 1.0, "bead_radius_a_over_b": 0.1},
        "body": {
            "prism": {
                "n_prism": 3,
                "dz_over_b": 0.5,
                "radius_over_b": 0.5,
                "axis": "x",
            },
            "length_total_um": 2.0,
        },
        "fluid": {"viscosity_Pa_s": 1.0e-3},
        "motor": {"torque_Nm": -1.0, "reverse_n_flagella": 1},
    }


def test_dt_over_tau_input_is_rejected() -> None:
    cfg = _base_cfg()
    cfg["time"] = {"duration_s": 0.1, "dt_over_tau": 0.001}

    with pytest.raises(ValueError, match="time.dt_over_tau"):
        SimulationConfig.from_dict(cfg)


def test_validate_time_scaling_is_always_fixed_to_paper_dt_star() -> None:
    cfg = _base_cfg()
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-7}
    sim_cfg = SimulationConfig.from_dict(cfg)

    sim_cfg.validate_time_scaling()
    assert sim_cfg.tau_s == pytest.approx(1.0)
    assert sim_cfg.dt_s == pytest.approx(1.0e-3)
    assert sim_cfg.dt_star == pytest.approx(1.0e-3)
    assert sim_cfg.output_dt_s == pytest.approx(1.0e-7)


def test_validate_time_scaling_for_motor_off_is_always_fixed() -> None:
    cfg = _base_cfg()
    cfg["motor"]["torque_Nm"] = 0.0
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-7}
    sim_cfg = SimulationConfig.from_dict(cfg)

    sim_cfg.validate_time_scaling()
    assert sim_cfg.tau_s == pytest.approx(1.0)
    assert sim_cfg.dt_s == pytest.approx(1.0e-3)
    assert sim_cfg.dt_star == pytest.approx(1.0e-3)


def test_validate_time_scaling_for_explicit_torque_is_always_fixed() -> None:
    cfg = _base_cfg()
    cfg["motor"]["torque_Nm"] = 1.0e-18
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-7}
    sim_cfg = SimulationConfig.from_dict(cfg)

    sim_cfg.validate_time_scaling()
    assert sim_cfg.tau_s == pytest.approx(1.0)
    assert sim_cfg.dt_s == pytest.approx(1.0e-3)
    assert sim_cfg.dt_star == pytest.approx(1.0e-3)


def test_torque_minus_one_uses_eta_b3_and_tau_is_one() -> None:
    cfg = _base_cfg()
    cfg["fluid"]["viscosity_Pa_s"] = 2.0e-3
    cfg["motor"]["torque_Nm"] = -1.0
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert sim_cfg.use_eta_b3_torque
    assert not sim_cfg.is_motor_off_torque
    assert sim_cfg.input_torque_Nm == pytest.approx(-1.0)
    assert sim_cfg.torque_scale_Nm == pytest.approx(
        sim_cfg.viscosity_Pa_s * (sim_cfg.b_m**3)
    )
    assert sim_cfg.torque_for_forces_Nm == pytest.approx(
        sim_cfg.viscosity_Pa_s * (sim_cfg.b_m**3)
    )
    assert sim_cfg.motor_torque_Nm == pytest.approx(
        sim_cfg.viscosity_Pa_s * (sim_cfg.b_m**3)
    )
    assert sim_cfg.torque_Nm == pytest.approx(sim_cfg.motor_torque_Nm)
    assert sim_cfg.tau_s == pytest.approx(1.0)


def test_torque_non_minus_one_uses_input_value() -> None:
    cfg = _base_cfg()
    cfg["motor"]["torque_Nm"] = 9.9e-18
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert not sim_cfg.use_eta_b3_torque
    assert not sim_cfg.is_motor_off_torque
    assert sim_cfg.input_torque_Nm == pytest.approx(9.9e-18)
    assert sim_cfg.torque_scale_Nm == pytest.approx(abs(9.9e-18))
    assert sim_cfg.torque_for_forces_Nm == pytest.approx(abs(9.9e-18))
    assert sim_cfg.motor_torque_Nm == pytest.approx(9.9e-18)
    assert sim_cfg.torque_Nm == pytest.approx(9.9e-18)
    assert sim_cfg.tau_s == pytest.approx(1.0)


def test_torque_for_forces_override_decouples_from_motor_torque() -> None:
    cfg = _base_cfg()
    cfg["motor"]["torque_Nm"] = 4.0e-21
    cfg["motor"]["torque_for_forces_override_Nm"] = 1.0e-21
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert sim_cfg.motor_torque_Nm == pytest.approx(4.0e-21)
    assert sim_cfg.torque_for_forces_Nm == pytest.approx(1.0e-21)


def test_torque_for_forces_override_non_positive_uses_default_coupling() -> None:
    cfg = _base_cfg()
    cfg["motor"]["torque_Nm"] = 4.0e-21
    cfg["motor"]["torque_for_forces_override_Nm"] = 0.0
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert sim_cfg.motor_torque_Nm == pytest.approx(4.0e-21)
    assert sim_cfg.torque_for_forces_Nm == pytest.approx(4.0e-21)


def test_torque_zero_sets_motor_off_but_keeps_tau_unity_scale() -> None:
    cfg = _base_cfg()
    cfg["motor"]["torque_Nm"] = 0.0
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert not sim_cfg.use_eta_b3_torque
    assert sim_cfg.is_motor_off_torque
    assert sim_cfg.motor_torque_Nm == pytest.approx(0.0)
    assert sim_cfg.torque_scale_Nm == pytest.approx(
        sim_cfg.viscosity_Pa_s * (sim_cfg.b_m**3)
    )
    assert sim_cfg.torque_for_forces_Nm == pytest.approx(
        sim_cfg.viscosity_Pa_s * (sim_cfg.b_m**3)
    )
    assert sim_cfg.tau_s == pytest.approx(1.0)
    assert sim_cfg.dt_s == pytest.approx(1.0e-3)
    assert sim_cfg.dt_star == pytest.approx(1.0e-3)


def test_motor_enable_switching_defaults_to_false() -> None:
    cfg = _base_cfg()
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert sim_cfg.motor.enable_switching is False


def test_motor_enable_switching_can_be_enabled() -> None:
    cfg = _base_cfg()
    cfg["motor"]["enable_switching"] = True
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert sim_cfg.motor.enable_switching is True


def test_motor_torque_ramp_defaults_to_disabled() -> None:
    cfg = _base_cfg()
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert sim_cfg.motor.torque_ramp_enabled is False
    assert sim_cfg.motor.torque_ramp_duration_s == pytest.approx(0.0)


def test_motor_torque_ramp_can_be_configured() -> None:
    cfg = _base_cfg()
    cfg["motor"]["torque_ramp_enabled"] = True
    cfg["motor"]["torque_ramp_duration_s"] = 0.05
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert sim_cfg.motor.torque_ramp_enabled is True
    assert sim_cfg.motor.torque_ramp_duration_s == pytest.approx(0.05)


def test_motor_local_scales_can_be_configured() -> None:
    cfg = _base_cfg()
    cfg["motor"]["local_hook_scale"] = 1.5
    cfg["motor"]["local_spring_scale"] = 1.25
    cfg["motor"]["local_bend_scale"] = 0.75
    cfg["motor"]["local_torsion_scale"] = 0.5
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert sim_cfg.motor.local_hook_scale == pytest.approx(1.5)
    assert sim_cfg.motor.local_spring_scale == pytest.approx(1.25)
    assert sim_cfg.motor.local_bend_scale == pytest.approx(0.75)
    assert sim_cfg.motor.local_torsion_scale == pytest.approx(0.5)


def test_projection_config_is_ignored_after_removal() -> None:
    cfg = _base_cfg()
    cfg["projection"] = {
        "enable_body_rigid_projection": True,
        "enable_hook_length_projection": True,
        "enable_basal_link_direction_projection": True,
        "enable_flagella_chain_length_projection": True,
        "enable_flagella_template_projection": True,
    }
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert not hasattr(sim_cfg, "projection")


def test_body_equiv_load_defaults_to_disabled() -> None:
    cfg = _base_cfg()
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert sim_cfg.body_equiv_load.enabled is False
    assert sim_cfg.body_equiv_load.mode == "none"
    assert sim_cfg.body_equiv_load.target_torque_Nm == pytest.approx(0.0)
    assert sim_cfg.body_equiv_load.target_force_N == pytest.approx(0.0)
    assert sim_cfg.body_equiv_load.attach_region_id == 0


def test_body_equiv_load_can_be_overridden() -> None:
    cfg = _base_cfg()
    cfg["body_equiv_load"] = {
        "enabled": True,
        "mode": "pure_couple",
        "target_torque_Nm": 5.0e-20,
        "target_force_N": 2.0e-13,
        "attach_region_id": 2,
    }
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert sim_cfg.body_equiv_load.enabled is True
    assert sim_cfg.body_equiv_load.mode == "pure_couple"
    assert sim_cfg.body_equiv_load.target_torque_Nm == pytest.approx(5.0e-20)
    assert sim_cfg.body_equiv_load.target_force_N == pytest.approx(2.0e-13)
    assert sim_cfg.body_equiv_load.attach_region_id == 2


def test_body_n_layers_is_derived_from_length_and_spacing() -> None:
    cfg = _base_cfg()
    cfg["body"]["length_total_um"] = 2.5
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)
    assert sim_cfg.compute_body_n_layers() == 6


def test_flagella_init_mode_defaults_and_bead_count_is_derived() -> None:
    cfg = _base_cfg()
    cfg["flagella"] = {"bond_L_over_b": 0.58, "length_over_b": 5.8}
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert sim_cfg.flagella.init_mode == "legacy_radius_pitch"
    assert sim_cfg.flagella.n_beads_per_flagellum == 11


def test_flagella_init_mode_can_be_set_to_paper_table1() -> None:
    cfg = _base_cfg()
    cfg["flagella"] = {
        "init_mode": "paper_table1",
        "n_beads_per_flagellum": 15,
        "bond_L_over_b": 0.58,
        "length_over_b": 5.8,
    }
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert sim_cfg.flagella.init_mode == "paper_table1"
    assert sim_cfg.flagella.n_beads_per_flagellum == 15


def test_flagella_stub_mode_can_be_set_to_minimal_basal_stub() -> None:
    cfg = _base_cfg()
    cfg["flagella"] = {
        "stub_mode": "minimal_basal_stub",
        "bond_L_over_b": 0.58,
        "length_over_b": 5.8,
    }
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert sim_cfg.flagella.stub_mode == "minimal_basal_stub"


def test_flagella_stub_mode_can_be_set_to_extended_basal_stub_5() -> None:
    cfg = _base_cfg()
    cfg["flagella"] = {
        "stub_mode": "extended_basal_stub_5",
        "bond_L_over_b": 0.58,
        "length_over_b": 5.8,
    }
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert sim_cfg.flagella.stub_mode == "extended_basal_stub_5"


def test_torsion_fd_eps_over_b_defaults_to_point_one() -> None:
    cfg = _base_cfg()
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert sim_cfg.potentials.torsion.fd_eps_over_b == pytest.approx(0.1)


def test_torsion_fd_eps_over_b_can_be_overridden() -> None:
    cfg = _base_cfg()
    cfg["potentials"] = {"torsion": {"fd_eps_over_b": 0.001}}
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert sim_cfg.potentials.torsion.fd_eps_over_b == pytest.approx(0.001)


def test_render2d_flagella_default_is_off() -> None:
    cfg = _base_cfg()
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert sim_cfg.render.render_flagella_2d is False


def test_render_center_body_in_2d_defaults_to_true() -> None:
    cfg = _base_cfg()
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert sim_cfg.render.center_body_in_2d is True


def test_render_center_body_in_2d_can_be_disabled() -> None:
    cfg = _base_cfg()
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    cfg["render"] = {"center_body_in_2d": False}
    sim_cfg = SimulationConfig.from_dict(cfg)

    assert sim_cfg.render.center_body_in_2d is False


def test_body_n_layers_requires_integer_multiple() -> None:
    cfg = _base_cfg()
    cfg["body"]["length_total_um"] = 2.3
    cfg["time"] = {"duration_s": 0.1, "dt_s": 1.0e-3}
    sim_cfg = SimulationConfig.from_dict(cfg)

    with pytest.raises(ValueError, match="整数倍"):
        sim_cfg.compute_body_n_layers()
