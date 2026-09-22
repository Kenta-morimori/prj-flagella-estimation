from __future__ import annotations

import json
from pathlib import Path

import pytest

from sim_swim.analysis.sweeps.generic_multi_run import run_campaign


class _ImmediateStop:
    """Test double: request a cooperative stop after the first internal step."""

    instances: list["_ImmediateStop"] = []

    def __init__(self) -> None:
        self.restored = False
        self.instances.append(self)

    def requested(self) -> str:
        return "interrupted by test SIGINT"

    def restore(self) -> None:
        self.restored = True


def test_generic_compact_interrupt_keeps_atomic_partial_evidence(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _ImmediateStop.instances.clear()
    campaign_path = tmp_path / "campaign.yaml"
    output_dir = tmp_path / "campaign-output"
    campaign_path.write_text(
        "\n".join(
            [
                "kind: generic_multi_run",
                "base_config: conf/sim_swim_2010.yaml",
                "base_overrides:",
                "  time.duration_s: 0.0001",
                "  output.policy: compact",
                "  output.checkpoint_interval_steps: 2500",
                "  motor.reference_torque_Nm: 2.0e-20",
                "sweep:",
                "  axes:",
                "    torque:",
                "      key: motor.torque_Nm",
                "      values: [2.0e-20]",
                "output:",
                f"  base_dir: {output_dir}",
                "  timestamp_subdir: false",
                "  save_state_archive: true",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    monkeypatch.setattr("sim_swim.core.run_context._require_clean_git", lambda: None)
    monkeypatch.setattr(
        "sim_swim.analysis.sweeps.generic_multi_run._StopRequest", _ImmediateStop
    )

    with pytest.raises(SystemExit) as exc_info:
        run_campaign(["--campaign-config", str(campaign_path)])
    assert exc_info.value.code == 130

    condition_dir = next(
        path
        for path in output_dir.iterdir()
        if path.is_dir() and (path / "progress.json").is_file()
    )
    progress = json.loads((condition_dir / "progress.json").read_text())
    summary = json.loads((condition_dir / "run_summary.json").read_text())
    performance = json.loads((condition_dir / "performance.json").read_text())
    completion = json.loads((output_dir / "campaign_completion.json").read_text())
    assert progress["status"] == "partial"
    assert progress["execution"]["completed_steps"] == 1
    samples_path = condition_dir / "diagnostic_samples.csv"
    assert samples_path.is_file()
    sample_header = samples_path.read_text(encoding="utf-8").splitlines()[0]
    assert "diagnostic_step" in sample_header
    assert "body_step" in sample_header
    assert (condition_dir / "state_archive.partial.npz").is_file()
    assert (condition_dir / "trajectory.partial.csv").is_file()
    assert summary["execution"]["status"] == "partial"
    assert performance["completed_steps"] == 1
    assert completion["status"] == "partial"
    assert completion["exit_code"] == 130
    assert not (output_dir / "summary.csv").exists()
    assert not (output_dir / "run_manifest.json").exists()
    assert _ImmediateStop.instances[0].restored is True


def test_generic_debug_campaign_does_not_install_cooperative_signal_handler(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    campaign_path = tmp_path / "debug-campaign.yaml"
    output_dir = tmp_path / "debug-output"
    campaign_path.write_text(
        "\n".join(
            [
                "kind: generic_multi_run",
                "base_config: conf/sim_swim_2010.yaml",
                "base_overrides:",
                "  time.duration_s: 5.0e-6",
                "  output.policy: debug",
                "  motor.reference_torque_Nm: 2.0e-20",
                "sweep:",
                "  axes:",
                "    torque:",
                "      key: motor.torque_Nm",
                "      values: [2.0e-20]",
                "output:",
                f"  base_dir: {output_dir}",
                "  timestamp_subdir: false",
                "  save_state_archive: false",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    class _ForbiddenStop:
        def __init__(self) -> None:
            raise AssertionError("debug campaign must not replace signal handlers")

    monkeypatch.setattr("sim_swim.core.run_context._require_clean_git", lambda: None)
    monkeypatch.setattr(
        "sim_swim.analysis.sweeps.generic_multi_run._StopRequest", _ForbiddenStop
    )

    run_campaign(["--campaign-config", str(campaign_path)])

    assert (output_dir / "summary.csv").is_file()
