#!/usr/bin/env python3
"""Persist and run a single, manifest-verified cs10 execution queue.

The dispatcher deliberately uses ``wait`` through :mod:`subprocess`; it never
polls Codex.  It is intended to be launched by a user in tmux on cs10.
"""

from __future__ import annotations

import argparse
import fcntl
import json
import os
from pathlib import Path
import shutil
import sqlite3
import subprocess
from typing import Callable
from datetime import UTC, datetime


DEFAULT_MIN_FREE_GIB = 10


def _now() -> str:
    return datetime.now(UTC).isoformat()


def _connect(state_dir: Path) -> sqlite3.Connection:
    state_dir.mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(state_dir / "queue.sqlite3")
    db.row_factory = sqlite3.Row
    db.execute(
        """CREATE TABLE IF NOT EXISTS jobs (
        id INTEGER PRIMARY KEY, branch TEXT NOT NULL, commit_sha TEXT NOT NULL,
        command TEXT NOT NULL, status TEXT NOT NULL, created_at TEXT NOT NULL,
        started_at TEXT, finished_at TEXT, worktree TEXT, output_root TEXT,
        detail TEXT)"""
    )
    return db


def _git(repo: Path, *args: str) -> str:
    return subprocess.run(
        ["git", *args], cwd=repo, check=True, text=True, capture_output=True
    ).stdout.strip()


def enqueue(db: sqlite3.Connection, repo: Path, branch: str, command: str) -> int:
    sha = _git(repo, "rev-parse", f"{branch}^{{commit}}")
    cursor = db.execute(
        "INSERT INTO jobs(branch, commit_sha, command, status, created_at) VALUES(?,?,?,?,?)",
        (branch, sha, command, "queued", _now()),
    )
    db.commit()
    return int(cursor.lastrowid)


def _manifest_succeeded(output_root: Path) -> tuple[bool, str]:
    path = output_root / "job_manifest.json"
    try:
        record = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return False, f"invalid or missing manifest: {exc}"
    if record.get("status") != "succeeded" or record.get("failed_configs"):
        return False, "job_manifest.json is not succeeded with empty failed_configs"
    return True, ""


def _notify(email: str | None, subject: str, body: str) -> None:
    if not email:
        return
    subprocess.run(
        ["/usr/bin/mail", "-s", subject, email], input=body, text=True, check=True
    )


def _free_gib(path: Path) -> float:
    usage = shutil.disk_usage(path)
    return usage.free / 1024**3


def run_once(
    db: sqlite3.Connection,
    repo: Path,
    state_dir: Path,
    email: str | None,
    min_free_gib: float = DEFAULT_MIN_FREE_GIB,
) -> int:
    job = db.execute(
        "SELECT * FROM jobs WHERE status='queued' ORDER BY id LIMIT 1"
    ).fetchone()
    if job is None:
        return 0
    free = _free_gib(state_dir)
    if free < min_free_gib:
        detail = f"free space {free:.2f} GiB is below {min_free_gib:.2f} GiB"
        db.execute(
            "UPDATE jobs SET status='blocked', detail=? WHERE id=?", (detail, job["id"])
        )
        db.commit()
        _notify(email, f"[cs10 queue] blocked job {job['id']}", detail)
        return 1
    worktree = state_dir / "worktrees" / f"job-{job['id']}-{job['commit_sha'][:12]}"
    worktree.parent.mkdir(parents=True, exist_ok=True)
    _git(repo, "worktree", "add", "--detach", str(worktree), job["commit_sha"])
    db.execute(
        "UPDATE jobs SET status='running', started_at=?, worktree=? WHERE id=?",
        (_now(), str(worktree), job["id"]),
    )
    db.commit()
    result = subprocess.run(
        job["command"], cwd=worktree, shell=True, text=True, capture_output=True
    )
    log = state_dir / f"job-{job['id']}.log"
    log.write_text(result.stdout + result.stderr, encoding="utf-8")
    output = next(
        (Path(line) for line in reversed(result.stdout.splitlines()) if line.strip()),
        None,
    )
    ok, detail = (False, f"command exited {result.returncode}")
    if result.returncode == 0 and output is not None:
        ok, detail = _manifest_succeeded(output)
    status = "succeeded" if ok else "failed"
    db.execute(
        "UPDATE jobs SET status=?, finished_at=?, output_root=?, detail=? WHERE id=?",
        (status, _now(), str(output) if output else None, detail, job["id"]),
    )
    db.commit()
    _notify(
        email,
        f"[cs10 queue] {status} job {job['id']}",
        f"branch={job['branch']}\ncommit={job['commit_sha']}\nlog={log}\n{detail}\n",
    )
    return 0 if ok else 1


def _with_dispatcher_lock(state_dir: Path, action: Callable[[], int]) -> int:
    """Run one dispatcher action while rejecting a second dispatcher."""
    lock_path = state_dir / "dispatcher.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("w", encoding="utf-8") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(
                "another cs10 queue dispatcher is already running"
            ) from exc
        return action()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument(
        "--state-dir",
        type=Path,
        default=Path(os.environ.get("XDG_STATE_HOME", Path.home() / ".local/state"))
        / "flagella-cs10-queue",
    )
    parser.add_argument("--email", default=os.environ.get("CS10_QUEUE_NOTIFY_EMAIL"))
    sub = parser.add_subparsers(dest="action", required=True)
    add = sub.add_parser("enqueue")
    add.add_argument("--branch", required=True)
    add.add_argument("--command", required=True)
    run = sub.add_parser("run-once")
    run.add_argument("--min-free-gib", type=float, default=DEFAULT_MIN_FREE_GIB)
    run_all = sub.add_parser("run-all")
    run_all.add_argument("--min-free-gib", type=float, default=DEFAULT_MIN_FREE_GIB)
    sub.add_parser("status")
    args = parser.parse_args(argv)
    db = _connect(args.state_dir)
    if args.action == "enqueue":
        print(enqueue(db, args.repo.resolve(), args.branch, args.command))
        return 0
    if args.action == "status":
        print(
            json.dumps(
                [dict(row) for row in db.execute("SELECT * FROM jobs ORDER BY id")],
                indent=2,
            )
        )
        return 0
    if not args.email:
        parser.error("CS10_QUEUE_NOTIFY_EMAIL or --email is required")
    if args.action == "run-once":
        return _with_dispatcher_lock(
            args.state_dir,
            lambda: run_once(
                db, args.repo.resolve(), args.state_dir, args.email, args.min_free_gib
            ),
        )

    def run_all() -> int:
        while db.execute("SELECT 1 FROM jobs WHERE status='queued' LIMIT 1").fetchone():
            result = run_once(
                db, args.repo.resolve(), args.state_dir, args.email, args.min_free_gib
            )
            if result:
                return result
        return 0

    return _with_dispatcher_lock(args.state_dir, run_all)


if __name__ == "__main__":
    raise SystemExit(main())
