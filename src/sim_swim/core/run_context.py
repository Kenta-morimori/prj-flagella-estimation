"""実行コンテキストの初期化ユーティリティ。"""

from __future__ import annotations

import json
import logging
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any
from zoneinfo import ZoneInfo


@dataclass(frozen=True)
class GitInfo:
    """実行時のGit情報。"""

    commit: str
    commit_short: str
    branch: str
    is_clean: bool


@dataclass(frozen=True)
class OutputPaths:
    """出力ディレクトリ群。"""

    root: Path
    sim_dir: Path
    render_dir: Path
    render2d_dir: Path


@dataclass(frozen=True)
class RunContext:
    """実行単位の共通情報。"""

    git: GitInfo
    out: OutputPaths
    logger: logging.Logger = field(repr=False)


def _run(cmd: list[str]) -> str:
    return subprocess.check_output(cmd, stderr=subprocess.STDOUT, text=True).strip()


def _git_info() -> GitInfo:
    """Git情報を取得する。Git管理外でも例外にせず unknown を返す。"""

    try:
        inside = _run(["git", "rev-parse", "--is-inside-work-tree"])
        if inside != "true":
            raise RuntimeError("not in git work tree")
        commit = _run(["git", "rev-parse", "HEAD"])
        commit_short = _run(["git", "rev-parse", "--short", "HEAD"])
        branch = _run(["git", "rev-parse", "--abbrev-ref", "HEAD"])
        status = _run(["git", "status", "--porcelain"])
        return GitInfo(
            commit=commit,
            commit_short=commit_short,
            branch=branch,
            is_clean=(status == ""),
        )
    except Exception:
        return GitInfo(
            commit="unknown",
            commit_short="unknown",
            branch="unknown",
            is_clean=False,
        )


def _now_jst() -> tuple[str, str]:
    dt = datetime.now(ZoneInfo("Asia/Tokyo"))
    return dt.strftime("%Y-%m-%d"), dt.strftime("%H%M%S")


def _setup_logger(log_path: Path) -> logging.Logger:
    logger = logging.getLogger("sim_swim")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False

    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)

    return logger


def init_run(base_dir: str | Path, input_info: dict[str, Any]) -> RunContext:
    """実行ディレクトリ・ログ・manifestを作成して返す。

    Args:
        base_dir: 実行結果の保存先ルート。
        input_info: 実行入力情報（configパス、overrideなど）。

    Returns:
        初期化済みの実行コンテキスト。
    """

    date_str, time_str = _now_jst()
    root = Path(base_dir) / date_str / time_str
    sim_dir = root / "sim"
    render_dir = root / "render"
    render2d_dir = root / "render2d"

    sim_dir.mkdir(parents=True, exist_ok=False)
    render_dir.mkdir(parents=True, exist_ok=False)
    render2d_dir.mkdir(parents=True, exist_ok=False)

    out = OutputPaths(
        root=root,
        sim_dir=sim_dir,
        render_dir=render_dir,
        render2d_dir=render2d_dir,
    )

    log_path = root / "run.log"
    logger = _setup_logger(log_path)
    git = _git_info()

    logger.info("Run start")
    logger.info("git.commit=%s", git.commit)
    logger.info("git.branch=%s", git.branch)
    logger.info("git.is_clean=%s", git.is_clean)
    logger.info("output_dir=%s", root)

    manifest = {
        "git": {
            "commit": git.commit,
            "commit_short": git.commit_short,
            "branch": git.branch,
            "is_clean": git.is_clean,
        },
        "run_time": {
            "date": date_str,
            "time": time_str,
            "timezone": "Asia/Tokyo",
        },
        "input": input_info,
        "outputs": {
            "root": str(root),
            "sim_dir": str(sim_dir),
            "render_dir": str(render_dir),
            "render2d_dir": str(render2d_dir),
            "log": str(log_path),
        },
        "environment": {
            "python_version": sys.version.split()[0],
        },
    }
    (root / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    return RunContext(git=git, out=out, logger=logger)
