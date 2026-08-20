from __future__ import annotations

import json
from pathlib import Path

import yaml

from sim_swim.analysis.partial_generic_multi_run import export_completed_campaign


def _run_summary(*, completed: bool) -> dict[str, object]:
    return {
        "execution": {
            "status": "completed" if completed else "running",
            "observed_final_t_s": 0.0 if completed else 0.00004,
            "expected_final_step_summary_t_s": 0.0,
        },
        "gates": {"shape_nonbody": {}, "shape_body": {}},
        "all_step_metrics": {},
    }


def test_export_completed_campaign_excludes_incomplete_conditions(
    tmp_path: Path,
) -> None:
    run_dir = tmp_path / "source"
    campaign_path = tmp_path / "campaign.yaml"
    campaign_path.write_text(
        yaml.safe_dump(
            {
                "base_config": "conf/sim_swim_2010.yaml",
                "base_overrides": {"time": {"duration_s": 0.0}},
                "sweep": {
                    "axes": {
                        "n_flagella": {
                            "key": "flagella.n_flagella",
                            "values": [1, 2],
                            "ids": ["nf01", "nf02"],
                        }
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    complete = run_dir / "nf01"
    complete.mkdir(parents=True)
    (complete / "run_summary.json").write_text(
        json.dumps(_run_summary(completed=True)), encoding="utf-8"
    )
    (complete / "state_archive.npz").write_bytes(b"archive")
    (complete / "trajectory.csv").write_text("t\n0\n", encoding="utf-8")

    incomplete = run_dir / "nf02"
    incomplete.mkdir()
    (incomplete / "run_summary.json").write_text(
        json.dumps(_run_summary(completed=False)), encoding="utf-8"
    )

    output_dir = export_completed_campaign(
        campaign_config=campaign_path,
        run_dir=run_dir,
        output_dir=tmp_path / "partial",
        overwrite=False,
    )

    manifest = json.loads((output_dir / "run_manifest.json").read_text())
    assert manifest["partial"] is True
    assert manifest["condition_order"] == ["nf01"]
    assert manifest["excluded_conditions"] == [
        {"condition_id": "nf02", "reason": "not_completed"}
    ]
    assert (output_dir / "summary.csv").is_file()
    assert (output_dir / "run.log").is_file()
