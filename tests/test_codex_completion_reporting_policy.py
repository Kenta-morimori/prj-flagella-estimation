from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
POLICY_PATHS = (
    ROOT / "AGENTS.md",
    ROOT / "tools/codex/skills/flagella-issue-workflow/SKILL.md",
    ROOT / "tools/codex/skills/flagella-issue-workflow/references/completion-policy.md",
    ROOT / "docs/codex/codex_workflow.md",
)


def test_completion_reporting_policy_requires_a_source_linked_pr() -> None:
    for path in POLICY_PATHS:
        policy = path.read_text(encoding="utf-8")
        assert "初回完了報告" in policy, path
        assert "PR作成後" in policy, path
        assert "PR不要" in policy, path
        assert "PR作成が失敗" in policy, path

    agents = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    assert "source Issue" in agents
    assert "Commentary progress updates remain allowed" in agents
