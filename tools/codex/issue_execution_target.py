#!/usr/bin/env python3
"""Map a GitHub Issue Form execution target to its canonical label."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


TARGET_TO_LABEL = {
    "mac_only": "execution:mac",
    "cs10_user_run": "execution:cs10",
    "no_runtime": "execution:none",
    "triage_required": "execution:triage",
}
HEADER = "### Heavy/runtime execution target"


def execution_label_from_issue_body(body: str | None) -> str:
    """Return the safe triage label when the Issue Form value is absent or invalid."""

    lines = (body or "").splitlines()
    try:
        header_index = lines.index(HEADER)
    except ValueError:
        return TARGET_TO_LABEL["triage_required"]
    for line in lines[header_index + 1 :]:
        value = line.strip()
        if not value:
            continue
        if value.startswith("### "):
            break
        return TARGET_TO_LABEL.get(value, TARGET_TO_LABEL["triage_required"])
    return TARGET_TO_LABEL["triage_required"]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--event", type=Path, required=True)
    args = parser.parse_args(argv)
    event = json.loads(args.event.read_text(encoding="utf-8"))
    body = dict(event.get("issue", {}) or {}).get("body")
    print(f"label={execution_label_from_issue_body(body)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
