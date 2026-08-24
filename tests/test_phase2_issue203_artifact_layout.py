import json
from pathlib import Path

from sim_swim.analysis.issue203_artifact_layout import (
    canonical_artifact_name,
    finalize_grid_composites,
    migrate_reference,
)


def _write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(value), encoding="utf-8")


def _root(tmp_path: Path, profile: str) -> Path:
    root = tmp_path / profile
    root.mkdir()
    conditions = []
    for n in (1, 2):
        condition_id = f"as000__ps001__nf{n:02d}"
        (root / condition_id).mkdir()
        conditions.append(
            {
                "condition_id": condition_id,
                "output_dir": f"/stale/{condition_id}",
                "axis_values": {"n_flagella": n},
            }
        )
    manifest = {"conditions": conditions}
    _write_json(root / "run_manifest.json", manifest)
    _write_json(root / "manifest.json", manifest)

    source_name = "composite_h264" if profile == "uniform" else "issue203_composite"
    source = root / "analysis" / source_name
    source.mkdir(parents=True)
    for record in conditions:
        condition_id = record["condition_id"]
        (source / f"{condition_id}_composite.mp4").touch()
        if profile == "diffusive":
            _write_json(
                source / f"{condition_id}_composite_manifest.json",
                {"condition_id": condition_id, "movie": "old.mp4"},
            )
    if profile == "uniform":
        legacy = root / "analysis" / "composite"
        legacy.mkdir()
        for record in conditions:
            condition_id = record["condition_id"]
            _write_json(
                legacy / f"{condition_id}_composite_manifest.json",
                {"condition_id": condition_id, "movie": "old.mp4"},
            )
        _write_json(source / "transcode_manifest.json", {"kind": "stale"})
    return root


def test_canonical_artifact_name_reorders_without_changing_condition_id() -> None:
    assert canonical_artifact_name("as002__ps001__nf03") == "nf03_as002_ps001"


def test_migration_dry_run_does_not_change_paths(tmp_path: Path, monkeypatch) -> None:
    root = _root(tmp_path, "uniform")
    monkeypatch.setattr(
        "sim_swim.analysis.issue203_artifact_layout._probe_video",
        lambda _: {"codec": "h264", "pix_fmt": "yuv420p", "frame_count": 21},
    )

    result = migrate_reference(root, "uniform")

    assert not result.applied
    assert (root / "as000__ps001__nf01").is_dir()
    assert (root / "analysis" / "composite_h264").is_dir()


def test_migration_applies_canonical_paths_and_manifests(
    tmp_path: Path, monkeypatch
) -> None:
    root = _root(tmp_path, "uniform")
    monkeypatch.setattr(
        "sim_swim.analysis.issue203_artifact_layout._probe_video",
        lambda _: {"codec": "h264", "pix_fmt": "yuv420p", "frame_count": 21},
    )

    result = migrate_reference(root, "uniform", apply=True)

    assert result.applied
    assert (root / "nf01_as000_ps001").is_dir()
    assert not (root / "as000__ps001__nf01").exists()
    assert (
        root / "analysis" / "composite" / "nf01_as000_ps001_composite.mp4"
    ).is_file()
    assert not (root / "analysis" / "composite_h264").exists()
    run_manifest = json.loads((root / "run_manifest.json").read_text())
    assert run_manifest["conditions"][0]["condition_id"] == "as000__ps001__nf01"
    assert run_manifest["conditions"][0]["output_dir"] == "nf01_as000_ps001"
    composite = json.loads(
        (
            root / "analysis" / "composite" / "nf01_as000_ps001_composite_manifest.json"
        ).read_text()
    )
    assert composite["movie"].endswith("nf01_as000_ps001_composite.mp4")
    assert composite["selected_codec"] == "libx264"


def test_finalize_grid_composites_removes_individual_artifacts(
    tmp_path: Path, monkeypatch
) -> None:
    root = tmp_path / "reference"
    composite = root / "analysis" / "composite"
    composite.mkdir(parents=True)
    for n in range(1, 4):
        (composite / f"nf{n:02d}_composite_grid.mp4").touch()
        _write_json(
            composite / f"nf{n:02d}_composite_grid_manifest.json",
            {"n_flagella": n, "condition_ids": [f"id_{i}" for i in range(9)]},
        )
    (composite / "nf01_as000_ps000_composite.mp4").touch()
    _write_json(composite / "nf01_as000_ps000_composite_manifest.json", {})
    monkeypatch.setattr(
        "sim_swim.analysis.issue203_artifact_layout._probe_video",
        lambda _: {"codec": "h264", "pix_fmt": "yuv420p", "frame_count": 21},
    )

    assert finalize_grid_composites(root, "uniform", apply=True) == 3

    assert not (composite / "nf01_as000_ps000_composite.mp4").exists()
    assert json.loads((composite / "manifest.json").read_text())["grid_count"] == 3
