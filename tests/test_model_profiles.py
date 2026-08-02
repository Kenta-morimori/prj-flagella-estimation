from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from sim_swim.sim.params import SimulationConfig

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
    ("sim_swim_2015.yaml", None),
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
def test_only_supported_project_profile_defaults_to_posterior_helix_axis(
    filename: str,
    expected: float | None,
) -> None:
    raw = yaml.safe_load((ROOT / "conf" / filename).read_text(encoding="utf-8"))
    cfg = SimulationConfig.from_dict(raw)

    assert cfg.flagella.initial_helix_axis_from_rear_deg == expected


def test_legacy_default_config_path_is_removed() -> None:
    assert not (ROOT / "conf" / "sim_swim.yaml").exists()
