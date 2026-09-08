from __future__ import annotations

from pathlib import Path

import yaml

from sim_swim.analysis.parallel_job import (
    build_plan,
    load_parallel_job,
    resolve_execution,
)


ROOT = Path(__file__).parents[1]


def test_agent_skill_and_runbook_share_the_cs10_parallel_gate() -> None:
    agent = (ROOT / "AGENTS.md").read_text(encoding="utf-8")
    skill = (ROOT / "tools/codex/skills/flagella-issue-workflow/SKILL.md").read_text(
        encoding="utf-8"
    )
    runbook = (ROOT / "docs/codex/cs10_runbook.md").read_text(encoding="utf-8")

    for phrase in (
        "two or more independent conditions",
        "cs10_qualified",
        "dry-run plan",
        "Issue-comment approval",
    ):
        assert phrase in agent
    for phrase in (
        "独立conditionが2以上",
        "cs10_qualified",
        "dry-run plan",
        "IssueコメントURL",
    ):
        assert phrase in skill
    for phrase in (
        "独立conditionが2以上",
        "cs10_qualified",
        "dry-run",
        "IssueコメントURL",
    ):
        assert phrase in runbook


def test_issue_form_requires_cs10_execution_contract_fields() -> None:
    form = yaml.safe_load(
        (ROOT / ".github/ISSUE_TEMPLATE/general.yml").read_text(encoding="utf-8")
    )
    fields = {item.get("id"): item for item in form["body"]}
    for field in (
        "cs10_execution_mode",
        "parallel_job_config",
        "parallel_worker_plan",
        "serial_exception_approval_url",
    ):
        assert fields[field]["validations"]["required"] is True
    assert fields["cs10_execution_mode"]["attributes"]["options"] == [
        "parallel",
        "serial_exception",
        "N/A",
    ]


def test_cs10_parallel_dry_run_exposes_workers_and_isolated_outputs(
    tmp_path: Path,
) -> None:
    job = load_parallel_job(
        ROOT / "conf/phase2_parallel/example_stage_a_validation/job.yaml"
    )
    execution = resolve_execution(job, None)
    plan = build_plan(job, execution, tmp_path / "parallel")

    assert execution.worker_policy == "cs10_qualified"
    assert execution.max_workers == 2
    assert len(plan["configs"]) == 2
    assert len({item["output_dir"] for item in plan["configs"]}) == 2
