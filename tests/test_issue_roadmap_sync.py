from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import yaml


ROOT = Path(__file__).resolve().parents[1]


def _tool():
    path = ROOT / "tools/codex/issue_roadmap_sync.py"
    spec = importlib.util.spec_from_file_location("issue_roadmap_sync", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _issue(
    tool, *, category: str | None, planned: str | None = None, closed: bool = False
):
    sections = []
    if category is not None:
        sections.extend([tool.ROADMAP_CATEGORY_HEADER, "", category])
    if planned is not None:
        sections.extend([tool.PLANNED_TARGET_DATE_HEADER, "", planned])
    return {
        "number": 217,
        "body": "\n\n".join(sections),
        "created_at": "2026-08-27T15:30:00Z",
        "closed_at": "2026-08-28T15:30:00Z" if closed else None,
        "labels": [],
    }


def test_issue_form_categories_match_milestone_mapping() -> None:
    tool = _tool()
    form = yaml.safe_load(
        (ROOT / ".github/ISSUE_TEMPLATE/general.yml").read_text(encoding="utf-8")
    )
    category = next(
        item for item in form["body"] if item.get("id") == "roadmap_category"
    )
    assert category["validations"]["required"] is True
    assert category["attributes"]["options"] == list(tool.MILESTONE_BY_CATEGORY)
    target_date = next(
        item for item in form["body"] if item.get("id") == "planned_target_date"
    )
    assert target_date.get("validations", {}).get("required") is not True

    workflow = (ROOT / ".github/workflows/issue-roadmap-sync.yml").read_text(
        encoding="utf-8"
    )
    assert "workflow_dispatch:" in workflow
    assert "PROJECT_AUTOMATION_TOKEN" in workflow
    assert "--dry-run" in workflow
    assert (
        "issue-roadmap-sync-${{ github.event.issue.number || inputs.issue_number }}"
        in workflow
    )


def test_opened_plan_sets_jst_start_and_optional_target() -> None:
    tool = _tool()
    category = "Phase 2→4 — Dataset & Model Selection"
    plan = tool.build_sync_plan(
        _issue(tool, category=category, planned="2026-09-30"), "opened"
    )
    assert plan.mode == "sync"
    assert plan.milestone == category
    assert plan.start_date == "2026-08-28"
    assert plan.planned_target_date == "2026-09-30"
    assert plan.close_date is None


def test_closed_plan_uses_jst_close_date_only_as_fallback() -> None:
    tool = _tool()
    plan = tool.build_sync_plan(
        _issue(tool, category="Phase 2 — Simulation & Pseudo-microscopy", closed=True),
        "closed",
    )
    assert plan.planned_target_date is None
    assert plan.close_date == "2026-08-29"


def test_no_response_optional_date_is_treated_as_unspecified() -> None:
    tool = _tool()
    plan = tool.build_sync_plan(
        _issue(
            tool,
            category="Phase 3 — Cell Clips",
            planned=tool.ISSUE_FORM_NO_RESPONSE,
        ),
        "opened",
    )
    assert plan.mode == "sync"
    assert plan.planned_target_date is None


def test_invalid_metadata_is_triaged_and_reopen_requires_review() -> None:
    tool = _tool()
    invalid_category = tool.build_sync_plan(_issue(tool, category=None), "opened")
    invalid_date = tool.build_sync_plan(
        _issue(tool, category="Phase 4 — Estimation Model", planned="2026/09/30"),
        "edited",
    )
    reopened = tool.build_sync_plan(_issue(tool, category=None), "reopened")
    assert invalid_category.mode == "triage"
    assert invalid_date.mode == "triage"
    assert reopened.mode == "needs_review"


def test_replace_roadmap_labels_preserves_other_label_namespaces() -> None:
    tool = _tool()

    class RecordingClient:
        def __init__(self) -> None:
            self.calls: list[tuple[str, str, dict | None]] = []

        def rest(self, method: str, path: str, payload: dict | None = None):
            self.calls.append((method, path, payload))
            if method == "GET" and path.endswith("/issues/217"):
                return {
                    "labels": [
                        {"name": "execution:none"},
                        {"name": "bug"},
                        {"name": tool.ROADMAP_TRIAGE_LABEL},
                    ]
                }
            return {}

    client = RecordingClient()
    tool.replace_roadmap_labels(
        client,
        "owner",
        "repo",
        {"number": 217},
        tool.ROADMAP_REVIEW_LABEL,
    )

    assert (
        "DELETE",
        "/repos/owner/repo/issues/217/labels/roadmap%3Atriage",
        None,
    ) in client.calls
    assert (
        "POST",
        "/repos/owner/repo/issues/217/labels",
        {"labels": [tool.ROADMAP_REVIEW_LABEL]},
    ) in client.calls
    assert not any(method == "PUT" for method, _, _ in client.calls)


def test_event_dry_run_needs_no_project_token(
    tmp_path: Path, monkeypatch, capsys
) -> None:
    tool = _tool()
    event = tmp_path / "event.json"
    event.write_text(
        json.dumps(
            {
                "action": "closed",
                "issue": _issue(
                    tool,
                    category="Phase 3→4 — Clip Dataset & Learning",
                    closed=True,
                ),
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.delenv("PROJECT_AUTOMATION_TOKEN", raising=False)
    assert tool.main(["--event", str(event), "--repo", "example", "--dry-run"]) == 0
    assert '"close_date": "2026-08-29"' in capsys.readouterr().out
