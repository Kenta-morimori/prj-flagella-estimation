from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import yaml


ROOT = Path(__file__).resolve().parents[1]


def _tool():
    path = ROOT / "tools/codex/issue_execution_target.py"
    spec = importlib.util.spec_from_file_location("issue_execution_target", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_issue_body_target_maps_to_exactly_one_execution_label() -> None:
    tool = _tool()
    for target, label in tool.TARGET_TO_LABEL.items():
        body = f"## 概要\n\n{tool.HEADER}\n\n{target}\n\n### 次\ntext"
        assert tool.execution_label_from_issue_body(body) == label


def test_missing_or_invalid_target_falls_back_to_triage() -> None:
    tool = _tool()
    assert tool.execution_label_from_issue_body(None) == "execution:triage"
    assert (
        tool.execution_label_from_issue_body(f"{tool.HEADER}\n\ninvalid")
        == "execution:triage"
    )


def test_event_cli_writes_github_output_label(tmp_path: Path, capsys) -> None:
    tool = _tool()
    event = tmp_path / "event.json"
    event.write_text(
        json.dumps({"issue": {"body": f"{tool.HEADER}\n\ncs10_user_run"}}),
        encoding="utf-8",
    )

    assert tool.main(["--event", str(event)]) == 0
    assert capsys.readouterr().out == "label=execution:cs10\n"


def test_issue_form_and_sync_workflow_use_the_same_execution_vocabulary() -> None:
    form = yaml.safe_load(
        (ROOT / ".github/ISSUE_TEMPLATE/general.yml").read_text(encoding="utf-8")
    )
    target = next(item for item in form["body"] if item.get("id") == "execution_target")
    assert target["validations"]["required"] is True
    assert target["attributes"]["options"] == [
        "mac_only",
        "cs10_user_run",
        "no_runtime",
        "triage_required",
    ]
    assert form["labels"] == ["execution:triage"]
    workflow = (ROOT / ".github/workflows/issue-execution-target.yml").read_text(
        encoding="utf-8"
    )
    assert "tools/codex/issue_execution_target.py" in workflow
    assert 'name.startsWith("execution:")' in workflow
    assert "github.rest.issues.removeLabel" in workflow
    assert "github.rest.issues.addLabels" in workflow
    assert "github.rest.issues.setLabels" not in workflow
    assert "issue-execution-label-${{ github.event.issue.number }}" in workflow
