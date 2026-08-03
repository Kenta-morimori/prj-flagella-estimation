from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from sim_swim.sim.params import SimulationConfig
from sim_swim.sim.core import Simulator

pytestmark = pytest.mark.light

ROOT = Path(__file__).resolve().parents[1]

PROFILE_CASES = (
    (
        "sim_swim_2010.yaml",
        (2010, "project", "legacy_project", "supported", 15, 11, 3, 48),
    ),
    (
        "sim_swim_2010_paper.yaml",
        (2010, "paper", "coarse", "supported", 15, 15, 3, 60),
    ),
    (
        "sim_swim_2015.yaml",
        (2015, "project", "refined", "pending", 30, 30, 3, 120),
    ),
    (
        "sim_swim_2015_paper.yaml",
        (2015, "paper", "refined", "pending", 30, 30, 3, 120),
    ),
)

INITIAL_HELIX_AXIS_CASES = (
    ("sim_swim_2010.yaml", 0.0),
    ("sim_swim_2010_paper.yaml", None),
    ("sim_swim_2015.yaml", 0.0),
    ("sim_swim_2015_paper.yaml", None),
)


@pytest.mark.parametrize(("filename", "expected"), PROFILE_CASES)
def test_canonical_model_profiles_have_expected_identity_and_resolution(
    filename: str,
    expected: tuple[int, str, str, str, int, int, int, int],
) -> None:
    raw = yaml.safe_load((ROOT / "conf" / filename).read_text(encoding="utf-8"))
    cfg = SimulationConfig.from_dict(raw)

    assert cfg.model_profile is not None
    profile = cfg.model_profile
    assert (
        profile.year,
        profile.variant,
        profile.resolution,
        profile.implementation_status,
        profile.body_beads,
        profile.flagellum_beads_per_filament,
        profile.nominal_flagella_count,
        profile.nominal_total_beads,
    ) == expected
    assert cfg.flagella.n_beads_per_flagellum == profile.flagellum_beads_per_filament
    assert cfg.flagella.n_flagella == profile.nominal_flagella_count
    assert cfg.body.prism.n_prism * cfg.compute_body_n_layers() == profile.body_beads
    assert profile.nominal_total_beads == profile.body_beads + (
        profile.flagellum_beads_per_filament * profile.nominal_flagella_count
    )


def test_2010_paper_profile_uses_fourteen_bonds_at_paper_spacing() -> None:
    raw = yaml.safe_load(
        (ROOT / "conf" / "sim_swim_2010_paper.yaml").read_text(encoding="utf-8")
    )
    cfg = SimulationConfig.from_dict(raw)

    assert cfg.flagella.n_beads_per_flagellum == 15
    assert cfg.flagella.bond_L_over_b == pytest.approx(0.58)
    assert cfg.flagella.length_over_b == pytest.approx(14 * 0.58)


@pytest.mark.parametrize(("filename", "expected"), INITIAL_HELIX_AXIS_CASES)
def test_project_profiles_default_to_posterior_helix_axis(
    filename: str,
    expected: float | None,
) -> None:
    raw = yaml.safe_load((ROOT / "conf" / filename).read_text(encoding="utf-8"))
    cfg = SimulationConfig.from_dict(raw)

    assert cfg.flagella.initial_helix_axis_from_rear_deg == expected


@pytest.mark.parametrize("filename", ["sim_swim_2015.yaml", "sim_swim_2015_paper.yaml"])
def test_2015_profiles_use_review_render_defaults(filename: str) -> None:
    raw = yaml.safe_load((ROOT / "conf" / filename).read_text(encoding="utf-8"))
    cfg = SimulationConfig.from_dict(raw)

    assert cfg.render.image_size_px == 96
    assert cfg.render.pixel_size_um == pytest.approx(0.1)
    assert cfg.render.view_range_um == pytest.approx(7.0)
    assert cfg.render.timestamp_fmt == "t = {t:.3f} s (τ = {tau_s:.3e} s)"
    assert cfg.render.projection_mode_2d == "body_capsule_orthographic_v1"


@pytest.mark.parametrize(
    ("filename", "force_distribution", "provenance"),
    [
        (
            "sim_swim_2015.yaml",
            "root_torque_segment_couples",
            "project_implementation",
        ),
        (
            "sim_swim_2015_paper.yaml",
            "hook_coupled_body_reaction",
            "paper_inspired_approximation",
        ),
    ],
)
def test_2015_profiles_are_geometry_implemented_and_evaluation_ready(
    filename: str,
    force_distribution: str,
    provenance: str,
) -> None:
    raw = yaml.safe_load((ROOT / "conf" / filename).read_text(encoding="utf-8"))
    cfg = SimulationConfig.from_dict(raw)

    assert cfg.motor.force_distribution == force_distribution
    implementation = cfg.implementation_manifest()
    expected_dynamics = {
        "implementation_status": "implemented",
        "force_distribution": force_distribution,
        "provenance": provenance,
    }
    if filename == "sim_swim_2015_paper.yaml":
        expected_dynamics.update(
            {
                "motor_axis_model": "hook_basal_tangent",
                "body_reaction_model": (
                    "local_attach_neighborhood_zero_net_force_torque_couple"
                ),
                "body_reaction_fallback_model": (
                    "all_body_beads_zero_net_force_torque_couple"
                ),
            }
        )
    assert implementation["dynamics"] == expected_dynamics
    assert implementation["geometry"] == {"implementation_status": "implemented"}
    assert implementation["simulation"] == {
        "implementation_status": "evaluation_ready",
        "blocked_by": [168],
    }
    cfg.validate_execution_supported()
    assert all(ok for _name, ok, _actual, _expected in cfg.paper_reference_checks())
    paper_reference = implementation["paper_reference"]
    assert paper_reference["parameters"]["body.width_over_b"] == {
        "value": 1.0,
        "source": "paper_default",
        "comparison_override": 0.7,
    }
    assert paper_reference["parameters"]["hook.length_over_b"] == {
        "value": 0.25,
        "source": "implementation_assumption",
        "implementation_status": "implemented",
    }
    assert cfg.flagella.placement_mode == "seeded_center_layer"
    assert cfg.body.prism.diagonal_braces_enabled is False


def test_2015_body_width_0p7_override_is_parseable_for_evaluation() -> None:
    raw = yaml.safe_load(
        (ROOT / "conf" / "sim_swim_2015_paper.yaml").read_text(encoding="utf-8")
    )
    cfg = SimulationConfig.from_dict(raw).with_overrides(
        {"body": {"prism": {"radius_over_b": 0.35}}}
    )

    assert 2.0 * cfg.body.prism.radius_over_b == pytest.approx(0.7)
    assert cfg.implementation_manifest()["paper_reference"]["parameters"][
        "body.width_over_b"
    ]["value"] == pytest.approx(0.7)
    cfg.validate_execution_supported()


def test_2015_paper_manifest_records_actual_refined_geometry() -> None:
    raw = yaml.safe_load(
        (ROOT / "conf" / "sim_swim_2015_paper.yaml").read_text(encoding="utf-8")
    )
    simulator = Simulator(SimulationConfig.from_dict(raw))
    manifest = simulator.implementation_manifest()
    actual = manifest["geometry"]["actual"]

    assert actual["total_beads"] == 120
    assert actual["body_beads"] == 30
    assert actual["flagellum_beads"] == [30, 30, 30]
    assert actual["body_layers"] == 5
    assert actual["body_slots_per_layer"] == 6
    assert actual["body_length_over_b"] == pytest.approx(2.0)
    assert actual["body_width_over_b"] == pytest.approx(1.0)
    assert actual["body_ring_edges"] == 30
    assert actual["body_vertical_edges"] == 24
    assert actual["body_diagonal_edges"] == 0
    assert actual["attachment_topology"] == [
        {"flag_id": 0, "body_bead_index": 12, "layer": 2, "slot": 0},
        {"flag_id": 1, "body_bead_index": 14, "layer": 2, "slot": 2},
        {"flag_id": 2, "body_bead_index": 16, "layer": 2, "slot": 4},
    ]


def test_non_2015_pending_profile_remains_blocked_in_manifest_and_execution() -> None:
    raw = yaml.safe_load(
        (ROOT / "conf" / "sim_swim_2010.yaml").read_text(encoding="utf-8")
    )
    raw["model_profile"]["implementation_status"] = "pending"
    cfg = SimulationConfig.from_dict(raw)

    assert cfg.implementation_manifest()["simulation"] == {
        "implementation_status": "blocked",
        "blocked_by": [],
    }
    with pytest.raises(ValueError, match="pending model profile"):
        cfg.validate_execution_supported()


def test_legacy_default_config_path_is_removed() -> None:
    assert not (ROOT / "conf" / "sim_swim.yaml").exists()
