from __future__ import annotations

from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]


def test_notification_workflow_requires_expected_inputs_and_gmail_secrets() -> None:
    workflow = ROOT / ".github/workflows/cs10-queue-notify.yml"
    document = yaml.load(workflow.read_text(encoding="utf-8"), Loader=yaml.BaseLoader)

    dispatch = document["on"]["workflow_dispatch"]
    inputs = dispatch["inputs"]
    assert inputs["event_type"]["required"] == "true"
    assert inputs["event_type"]["options"] == [
        "succeeded",
        "failed",
        "cancelled",
        "queue_completed",
    ]
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
    for secret in (
        "CS10_QUEUE_NOTIFY_EMAIL",
        "CS10_QUEUE_GMAIL_USERNAME",
        "CS10_QUEUE_GMAIL_APP_PASSWORD",
    ):
        assert f"secrets.{secret}" in source
    assert "smtp.gmail.com" in source
