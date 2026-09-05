from __future__ import annotations

from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]


def test_notification_workflow_records_events_without_mail_secrets() -> None:
    workflow = ROOT / ".github/workflows/cs10-queue-notify.yml"
    document = yaml.load(workflow.read_text(encoding="utf-8"), Loader=yaml.BaseLoader)

    dispatch = document["on"]["workflow_dispatch"]
    inputs = dispatch["inputs"]
    assert inputs["event_type"]["required"] == "true"
    assert inputs["event_type"]["options"] == ["job_completed"]
    assert set(inputs) == {
        "event_type",
        "reservation_id",
        "branch",
        "commit_sha",
        "state",
        "logs",
        "control_dir",
        "output_root",
        "error",
    }

    source = workflow.read_text(encoding="utf-8")
    assert "GITHUB_STEP_SUMMARY" in source
    assert 'if [ "$STATE" = "failed" ]; then' in source
    assert "exit 1" in source
    assert "secrets." not in source
    assert "smtp.gmail.com" not in source
    assert "GMAIL_" not in source
