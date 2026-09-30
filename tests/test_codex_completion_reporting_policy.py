from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_completion_reporting_policy_requires_a_source_linked_pr() -> None:
    for path in (
        ROOT / "AGENTS.md",
        ROOT / "tools/codex/skills/flagella-issue-workflow/SKILL.md",
        ROOT / "docs/codex/codex_workflow.md",
    ):
        policy = path.read_text(encoding="utf-8")
        assert "初回完了報告" in policy, path
        assert "source Issue" in policy, path

    agents = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    assert "PR 不要" in agents
    assert "PR 作成失敗" in agents
    assert "This file defines" not in agents

    completion_reference = (
        ROOT
        / "tools/codex/skills/flagella-issue-workflow/references/completion-policy.md"
    ).read_text(encoding="utf-8")
    assert "PR 不要" in completion_reference
    assert "PR 作成失敗" in completion_reference


def test_post_merge_handoff_preserves_parent_and_requires_approval_for_followups() -> (
    None
):
    workflow = (ROOT / "docs/codex/codex_workflow.md").read_text(encoding="utf-8")

    for phrase in (
        "PRマージ後のIssue引継ぎ",
        "source Issue",
        "受入条件",
        "review_result.json",
        "継続親Issue",
        "ユーザー承認後",
        "Heavy/runtime execution target",
    ):
        assert phrase in workflow
