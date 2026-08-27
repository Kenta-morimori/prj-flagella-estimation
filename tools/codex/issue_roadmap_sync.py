#!/usr/bin/env python3
"""Synchronize Issue Form roadmap metadata to the repository Project.

The module keeps parsing and date policy independent from GitHub API calls so the
workflow policy is unit-testable without credentials.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from datetime import date, datetime
import json
import os
from pathlib import Path
import re
from typing import Any
from urllib.error import HTTPError
from urllib.request import Request, urlopen
from zoneinfo import ZoneInfo


PROJECT_OWNER = "Kenta-morimori"
PROJECT_NUMBER = 8
JST = ZoneInfo("Asia/Tokyo")
ROADMAP_CATEGORY_HEADER = "### Roadmap category (Milestone)"
PLANNED_TARGET_DATE_HEADER = "### Planned target date"

MILESTONE_BY_CATEGORY = {
    "Phase 1 — Foundations": "Phase 1 — Foundations",
    "Phase 2 — Simulation & Pseudo-microscopy": "Phase 2 — Simulation & Pseudo-microscopy",
    "Phase 2→3 — Video Handoff": "Phase 2→3 — Video Handoff",
    "Phase 2→4 — Dataset & Model Selection": "Phase 2→4 — Dataset & Model Selection",
    "Phase 3 — Cell Clips": "Phase 3 — Cell Clips",
    "Phase 3→4 — Clip Dataset & Learning": "Phase 3→4 — Clip Dataset & Learning",
    "Phase 4 — Estimation Model": "Phase 4 — Estimation Model",
    "Project & Operations": "Project & Operations",
}
ROADMAP_TRIAGE_LABEL = "roadmap:triage"
ROADMAP_REVIEW_LABEL = "roadmap:needs-review"


@dataclass(frozen=True)
class SyncPlan:
    mode: str
    milestone: str | None = None
    start_date: str | None = None
    planned_target_date: str | None = None
    close_date: str | None = None
    triage_reason: str | None = None


def form_value(body: str | None, header: str) -> str | None:
    """Return the first non-empty Issue Form value below ``header``."""

    lines = (body or "").splitlines()
    try:
        header_index = lines.index(header)
    except ValueError:
        return None
    for line in lines[header_index + 1 :]:
        value = line.strip()
        if not value:
            continue
        if value.startswith("### "):
            return None
        return value
    return None


def jst_date(iso_timestamp: str) -> str:
    timestamp = datetime.fromisoformat(iso_timestamp.replace("Z", "+00:00"))
    return timestamp.astimezone(JST).date().isoformat()


def valid_date(value: str) -> str | None:
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", value):
        return None
    try:
        return date.fromisoformat(value).isoformat()
    except ValueError:
        return None


def build_sync_plan(issue: dict[str, Any], event_action: str) -> SyncPlan:
    """Build the non-mutating synchronization plan for an issue event."""

    if event_action == "reopened":
        return SyncPlan(mode="needs_review")

    category = form_value(issue.get("body"), ROADMAP_CATEGORY_HEADER)
    milestone = MILESTONE_BY_CATEGORY.get(category or "")
    if milestone is None:
        return SyncPlan(
            mode="triage", triage_reason="missing or invalid roadmap category"
        )

    planned_target = form_value(issue.get("body"), PLANNED_TARGET_DATE_HEADER)
    if planned_target and (normalized_target := valid_date(planned_target)) is None:
        return SyncPlan(mode="triage", triage_reason="invalid planned target date")

    start_date = jst_date(issue["created_at"])
    if event_action == "closed":
        closed_at = issue.get("closed_at")
        if not closed_at:
            return SyncPlan(
                mode="triage", triage_reason="closed issue has no closed_at"
            )
        return SyncPlan(
            mode="sync",
            milestone=milestone,
            start_date=start_date,
            planned_target_date=normalized_target if planned_target else None,
            close_date=jst_date(closed_at),
        )
    return SyncPlan(
        mode="sync",
        milestone=milestone,
        start_date=start_date,
        planned_target_date=normalized_target if planned_target else None,
    )


class GitHubClient:
    def __init__(self, token: str) -> None:
        self.token = token

    def rest(
        self, method: str, path: str, payload: dict[str, Any] | None = None
    ) -> Any:
        data = json.dumps(payload).encode() if payload is not None else None
        request = Request(
            f"https://api.github.com{path}",
            data=data,
            method=method,
            headers={
                "Accept": "application/vnd.github+json",
                "Authorization": f"Bearer {self.token}",
                "X-GitHub-Api-Version": "2022-11-28",
                "Content-Type": "application/json",
            },
        )
        with urlopen(request) as response:  # noqa: S310 - fixed GitHub API origin
            return json.load(response) if response.readable() else None

    def graphql(self, query: str, variables: dict[str, Any]) -> dict[str, Any]:
        result = self.rest("POST", "/graphql", {"query": query, "variables": variables})
        if result.get("errors"):
            raise RuntimeError(json.dumps(result["errors"], ensure_ascii=False))
        return result["data"]


def ensure_label(
    client: GitHubClient, owner: str, repo: str, name: str, color: str
) -> None:
    try:
        client.rest("GET", f"/repos/{owner}/{repo}/labels/{name}")
    except HTTPError as error:
        if error.code != 404:
            raise
        client.rest(
            "POST",
            f"/repos/{owner}/{repo}/labels",
            {
                "name": name,
                "color": color,
                "description": "Roadmap metadata requires review",
            },
        )


def replace_roadmap_labels(
    client: GitHubClient, owner: str, repo: str, issue: dict[str, Any], add: str | None
) -> None:
    current_issue = client.rest(
        "GET", f"/repos/{owner}/{repo}/issues/{issue['number']}"
    )
    labels = [label["name"] for label in current_issue.get("labels", [])]
    retained = [
        label
        for label in labels
        if label not in {ROADMAP_TRIAGE_LABEL, ROADMAP_REVIEW_LABEL}
    ]
    if add:
        ensure_label(client, owner, repo, add, "fbca04")
        retained.append(add)
    client.rest(
        "PUT",
        f"/repos/{owner}/{repo}/issues/{issue['number']}/labels",
        {"labels": retained},
    )


PROJECT_QUERY = """
query($owner: String!, $number: Int!, $issueNumber: Int!, $repo: String!) {
  user(login: $owner) {
    projectV2(number: $number) {
      id
      fields(first: 50) { nodes { ... on ProjectV2Field { id name } } }
    }
  }
  repository(owner: $owner, name: $repo) {
    issue(number: $issueNumber) {
      id
      projectItems(first: 50) {
        nodes {
          id
          project { id }
          fieldValues(first: 50) {
            nodes {
              ... on ProjectV2ItemFieldDateValue { date field { ... on ProjectV2Field { name } } }
            }
          }
        }
      }
    }
  }
}
"""


def project_item_and_fields(
    client: GitHubClient, owner: str, repo: str, issue_number: int
) -> tuple[str, str, str, str | None, str | None]:
    data = client.graphql(
        PROJECT_QUERY,
        {
            "owner": owner,
            "number": PROJECT_NUMBER,
            "repo": repo,
            "issueNumber": issue_number,
        },
    )
    project = data["user"]["projectV2"]
    if project is None:
        raise RuntimeError(f"Project #{PROJECT_NUMBER} is not accessible")
    fields = {
        field["name"]: field["id"] for field in project["fields"]["nodes"] if field
    }
    start_field = fields.get("Start date")
    target_field = fields.get("Target date")
    if not start_field or not target_field:
        raise RuntimeError("Project must provide Start date and Target date fields")

    issue = data["repository"]["issue"]
    items = [
        item
        for item in issue["projectItems"]["nodes"]
        if item["project"]["id"] == project["id"]
    ]
    if not items:
        added = client.graphql(
            "mutation($projectId: ID!, $contentId: ID!) { addProjectV2ItemById(input: {projectId: $projectId, contentId: $contentId}) { item { id } } }",
            {"projectId": project["id"], "contentId": issue["id"]},
        )
        return (
            project["id"],
            added["addProjectV2ItemById"]["item"]["id"],
            start_field,
            None,
            None,
        )

    item = items[0]
    dates = {
        value["field"]["name"]: value["date"]
        for value in item["fieldValues"]["nodes"]
        if value.get("field") and "date" in value
    }
    return (
        project["id"],
        item["id"],
        start_field,
        dates.get("Start date"),
        dates.get("Target date"),
    )


def set_project_date(
    client: GitHubClient, project_id: str, item_id: str, field_id: str, value: str
) -> None:
    client.graphql(
        "mutation($projectId: ID!, $itemId: ID!, $fieldId: ID!, $date: Date!) { updateProjectV2ItemFieldValue(input: {projectId: $projectId, itemId: $itemId, fieldId: $fieldId, value: {date: $date}}) { projectV2Item { id } } }",
        {
            "projectId": project_id,
            "itemId": item_id,
            "fieldId": field_id,
            "date": value,
        },
    )


def synchronize(
    client: GitHubClient, owner: str, repo: str, issue: dict[str, Any], action: str
) -> SyncPlan:
    plan = build_sync_plan(issue, action)
    if plan.mode == "triage":
        replace_roadmap_labels(client, owner, repo, issue, ROADMAP_TRIAGE_LABEL)
        return plan
    if plan.mode == "needs_review":
        replace_roadmap_labels(client, owner, repo, issue, ROADMAP_REVIEW_LABEL)
        return plan

    milestones = client.rest(
        "GET", f"/repos/{owner}/{repo}/milestones?state=all&per_page=100"
    )
    milestone = next(
        (entry for entry in milestones if entry["title"] == plan.milestone), None
    )
    if milestone is None:
        raise RuntimeError(f"Milestone not found: {plan.milestone}")
    client.rest(
        "PATCH",
        f"/repos/{owner}/{repo}/issues/{issue['number']}",
        {"milestone": milestone["number"]},
    )
    replace_roadmap_labels(client, owner, repo, issue, None)

    project_id, item_id, start_field, existing_start, existing_target = (
        project_item_and_fields(client, owner, repo, issue["number"])
    )
    if existing_start is None:
        set_project_date(
            client, project_id, item_id, start_field, plan.start_date or ""
        )

    fields = client.graphql(
        "query($projectId: ID!) { node(id: $projectId) { ... on ProjectV2 { fields(first: 50) { nodes { ... on ProjectV2Field { id name } } } } } }",
        {"projectId": project_id},
    )
    target_field = next(
        field["id"]
        for field in fields["node"]["fields"]["nodes"]
        if field and field["name"] == "Target date"
    )
    target = plan.planned_target_date or (
        plan.close_date if existing_target is None else None
    )
    if target:
        set_project_date(client, project_id, item_id, target_field, target)
    return plan


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--event", type=Path)
    parser.add_argument("--issue-number", type=int)
    parser.add_argument("--owner", default=PROJECT_OWNER)
    parser.add_argument("--repo", required=True)
    parser.add_argument("--action")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    if args.event:
        event = json.loads(args.event.read_text(encoding="utf-8"))
        issue = event["issue"]
        action = args.action or event["action"]
    else:
        if args.issue_number is None or args.action is None:
            raise SystemExit("--issue-number and --action are required without --event")
        token = os.environ.get("PROJECT_AUTOMATION_TOKEN")
        if not token:
            raise SystemExit(
                "PROJECT_AUTOMATION_TOKEN is required to load an issue by number"
            )
        client = GitHubClient(token)
        issue = client.rest(
            "GET", f"/repos/{args.owner}/{args.repo}/issues/{args.issue_number}"
        )
        action = args.action

    plan = build_sync_plan(issue, action)
    print(json.dumps(plan.__dict__, ensure_ascii=False, sort_keys=True))
    if args.dry_run:
        return 0
    token = os.environ.get("PROJECT_AUTOMATION_TOKEN")
    if not token:
        raise SystemExit(
            "PROJECT_AUTOMATION_TOKEN is required for Project #8 synchronization"
        )
    client = GitHubClient(token)
    synchronize(client, args.owner, args.repo, issue, action)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
