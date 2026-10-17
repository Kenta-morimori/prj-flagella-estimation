"""The common preview accepts campaigns independent of their mechanics labels."""

import importlib.util
from pathlib import Path

import yaml

from sim_swim.analysis.multi_run_campaign import load_yaml

ROOT = Path(__file__).resolve().parents[1]


def test_preview_full_vector_only_deduplicates_geometry_without_axis_name_presets(
    tmp_path,
):
    spec = importlib.util.spec_from_file_location(
        "initial_geometry_preview",
        ROOT / "scripts/01_simulate_swimming/preview_initial_geometry.py",
    )
    preview = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(preview)
    campaign = load_yaml(
        ROOT / "conf/phase2_multi_run/2010_hex_project_nf3_bundle_1tau_issue255.yaml"
    )
    axis = campaign["sweep"]["axes"]["attachment_pattern"]
    # Opaque IDs and count-key aliases must not control geometry selection.
    axis["keys"]["filaments"] = axis["keys"].pop("n_flagella")
    axis["keys"]["body_reaction_full_vector"] = "motor.body_reaction_full_vector"
    for value in axis["values"]:
        value["filaments"] = value.pop("n_flagella")
        value["body_reaction_full_vector"] = True
    axis["ids"] = ["alpha", "beta", "gamma", "delta"]
    config = tmp_path / "arbitrary_campaign.yaml"
    config.write_text(yaml.safe_dump(campaign))
    output = tmp_path / "preview"
    output.mkdir()
    report = preview.render_config(config, output)
    assert report["condition_count"] == 4
    assert report["unique_shape_count"] == 1
    assert report["overview_condition_ids"] == ["alpha"]
    assert len(report["spatial_limits_um"]) == 3
    for key in ("image", "axial_image", "overview_image"):
        assert Path(report[key]).is_file()
