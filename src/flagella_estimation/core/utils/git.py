from __future__ import annotations

import subprocess
from dataclasses import dataclass


@dataclass(frozen=True)
class GitInfo:
    commit: str  # full sha
    commit_short: str  # short sha
    branch: str


class GitNotReadyError(RuntimeError):
    """Raised when repository state is not suitable for a reproducible run."""


def _run(cmd: list[str]) -> str:
    return subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True).strip()


def require_committed_state() -> GitInfo:
    # Git管理下か
    try:
        inside = _run(["git", "rev-parse", "--is-inside-work-tree"])
    except Exception as e:
        raise GitNotReadyError("Git管理下ではない。git init してコミットせよ。") from e
    if inside != "true":
        raise GitNotReadyError("Git管理下ではない。git init してコミットせよ。")

    # 未コミット変更がないか
    status = _run(["git", "status", "--porcelain"])
    if status:
        raise GitNotReadyError(
            "未コミットの変更があるため実行不可である。commit（またはstash）してから実行せよ。"
        )

    commit = _run(["git", "rev-parse", "HEAD"])
    branch = _run(["git", "rev-parse", "--abbrev-ref", "HEAD"])
    commit_short = _run(["git", "rev-parse", "--short", "HEAD"])
    return GitInfo(commit=commit, commit_short=commit_short, branch=branch)
