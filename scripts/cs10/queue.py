#!/usr/bin/env python3
"""Persistent, single-worker cs10 queue for fixed-commit parallel jobs.

Run this program from a tmux session on cs10.  Reservations never execute an
arbitrary shell command: they resolve one repository branch to a commit and run
one checked-in parallel-job config through ``parallel_tmux.py run``.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime
import fcntl
import json
import os
from pathlib import Path
import signal
import sqlite3
import subprocess
import time
from typing import Any, Iterator
from uuid import uuid4
from zoneinfo import ZoneInfo


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_STATE_DIR = Path.home() / ".local/state/prj-flagella-estimation/cs10-queue"
DEFAULT_WORKTREE_DIR = Path.home() / "src/prj-flagella-estimation-queue-worktrees"
TERMINAL_STATES = {"succeeded", "failed", "cancelled", "blocked"}


def _now() -> str:
    return datetime.now(ZoneInfo("Asia/Tokyo")).isoformat()


def _state_dir() -> Path:
    return Path(os.environ.get("CS10_QUEUE_STATE_DIR", DEFAULT_STATE_DIR))


def _worktree_dir() -> Path:
    return Path(os.environ.get("CS10_QUEUE_WORKTREE_DIR", DEFAULT_WORKTREE_DIR))


@dataclass(frozen=True)
class Reservation:
    id: int
    branch: str
    commit_sha: str
    config: str
    priority: int
    state: str
    created_at: str
    worktree: str | None
    control_dir: str | None
    output_root: str | None
    pid: int | None
    exit_code: int | None
    error: str | None


class QueueStore:
    def __init__(self, state_dir: Path | None = None) -> None:
        self.state_dir = state_dir or _state_dir()
        self.state_dir.mkdir(parents=True, exist_ok=True)
        self.connection = sqlite3.connect(self.state_dir / "queue.sqlite3")
        self.connection.row_factory = sqlite3.Row
        self.connection.execute("PRAGMA foreign_keys=ON")
        self.connection.execute(
            """CREATE TABLE IF NOT EXISTS reservations (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                branch TEXT NOT NULL,
                commit_sha TEXT NOT NULL,
                config TEXT NOT NULL,
                priority INTEGER NOT NULL,
                state TEXT NOT NULL,
                created_at TEXT NOT NULL,
                worktree TEXT,
                control_dir TEXT,
                output_root TEXT,
                pid INTEGER,
                exit_code INTEGER,
                error TEXT
            )"""
        )
        self.connection.execute(
            "CREATE TABLE IF NOT EXISTS queue_control (id INTEGER PRIMARY KEY CHECK (id=1), paused INTEGER NOT NULL)"
        )
        self.connection.execute(
            "INSERT OR IGNORE INTO queue_control (id, paused) VALUES (1, 0)"
        )
        self.connection.execute(
            """CREATE TABLE IF NOT EXISTS events (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                reservation_id INTEGER,
                created_at TEXT NOT NULL,
                kind TEXT NOT NULL,
                detail TEXT NOT NULL,
                FOREIGN KEY(reservation_id) REFERENCES reservations(id)
            )"""
        )
        self.connection.commit()

    def close(self) -> None:
        self.connection.close()

    def event(self, reservation_id: int | None, kind: str, detail: str) -> None:
        self.connection.execute(
            "INSERT INTO events (reservation_id, created_at, kind, detail) VALUES (?, ?, ?, ?)",
            (reservation_id, _now(), kind, detail),
        )
        self.connection.commit()

    def add(
        self, *, branch: str, commit_sha: str, config: str, priority: int
    ) -> Reservation:
        cursor = self.connection.execute(
            "INSERT INTO reservations (branch, commit_sha, config, priority, state, created_at) VALUES (?, ?, ?, ?, 'queued', ?)",
            (branch, commit_sha, config, priority, _now()),
        )
        self.connection.commit()
        reservation = self.get(int(cursor.lastrowid))
        assert reservation is not None
        self.event(reservation.id, "queued", f"{branch}@{commit_sha} {config}")
        return reservation

    def get(self, reservation_id: int) -> Reservation | None:
        row = self.connection.execute(
            "SELECT * FROM reservations WHERE id=?", (reservation_id,)
        ).fetchone()
        return Reservation(**dict(row)) if row else None

    def list(self) -> list[Reservation]:
        rows = self.connection.execute(
            "SELECT * FROM reservations ORDER BY CASE state WHEN 'running' THEN 0 WHEN 'queued' THEN 1 ELSE 2 END, priority DESC, id ASC"
        ).fetchall()
        return [Reservation(**dict(row)) for row in rows]

    def next(self) -> Reservation | None:
        row = self.connection.execute(
            "SELECT * FROM reservations WHERE state='queued' ORDER BY priority DESC, id ASC LIMIT 1"
        ).fetchone()
        return Reservation(**dict(row)) if row else None

    def update(self, reservation_id: int, **values: Any) -> None:
        if not values:
            return
        columns = ", ".join(f"{column}=?" for column in values)
        self.connection.execute(
            f"UPDATE reservations SET {columns} WHERE id=?",
            (*values.values(), reservation_id),
        )
        self.connection.commit()

    def paused(self) -> bool:
        return bool(
            self.connection.execute(
                "SELECT paused FROM queue_control WHERE id=1"
            ).fetchone()[0]
        )

    def set_paused(self, paused: bool) -> None:
        self.connection.execute(
            "UPDATE queue_control SET paused=? WHERE id=1", (int(paused),)
        )
        self.connection.commit()
        self.event(
            None, "queue_paused" if paused else "queue_resumed", "operator action"
        )


def _git(*args: str, cwd: Path = REPOSITORY_ROOT) -> str:
    return subprocess.run(
        ["git", *args], cwd=cwd, check=True, text=True, capture_output=True
    ).stdout.strip()


def resolve_commit(branch: str) -> str:
    return _git("rev-parse", f"{branch}^{{commit}}")


def validate_config(config: str) -> str:
    path = Path(config)
    if path.is_absolute() or ".." in path.parts or path.suffix not in {".yaml", ".yml"}:
        raise ValueError("config must be a repository-relative YAML path")
    resolved = (REPOSITORY_ROOT / path).resolve()
    if not resolved.is_file() or REPOSITORY_ROOT not in resolved.parents:
        raise ValueError(f"parallel-job config is unavailable: {config}")
    return str(path)


def ensure_worktree(reservation: Reservation) -> Path:
    root = _worktree_dir()
    root.mkdir(parents=True, exist_ok=True)
    path = root / f"queue-{reservation.id:05d}-{reservation.commit_sha[:12]}"
    if path.exists():
        return path
    _git("worktree", "add", "--detach", str(path), reservation.commit_sha)
    return path


def _runtime_python() -> Path:
    python = REPOSITORY_ROOT / ".venv-cs10/bin/python"
    if not python.is_file():
        raise RuntimeError(f"missing shared cs10 runtime: {python}")
    return python


def _reservation_paths(
    store: QueueStore, reservation: Reservation
) -> tuple[Path, Path]:
    directory = store.state_dir / f"reservation-{reservation.id:05d}"
    directory.mkdir(exist_ok=True)
    return directory / "stdout.log", directory / "stderr.log"


def _load_launch_record(stdout_log: Path) -> dict[str, Any]:
    text = stdout_log.read_text(encoding="utf-8").strip()
    if not text:
        raise RuntimeError("foreground launcher did not write its launch record")
    return json.loads(text)


def _manifest_succeeded(record: dict[str, Any]) -> tuple[bool, str]:
    root = Path(str(record["output_root"]))
    manifest = root / "job_manifest.json"
    completion = root / "campaign/campaign_completion.json"
    if not manifest.is_file() or not completion.is_file():
        return False, "missing job manifest or campaign completion marker"
    job = json.loads(manifest.read_text(encoding="utf-8"))
    campaign = json.loads(completion.read_text(encoding="utf-8"))
    failed = job.get("failed_configs", [])
    records = list(job.get("configs", []) or [])
    expected = int(campaign.get("expected_condition_count", -1))
    conditions = root / "campaign/conditions"
    summaries = (
        sum((item / "run_summary.json").is_file() for item in conditions.iterdir())
        if conditions.is_dir()
        else 0
    )
    if (
        job.get("status") == "succeeded"
        and not failed
        and job.get("aggregation", {}).get("status") == "completed"
        and campaign.get("status") == "completed"
        and campaign.get("exit_code") == 0
        and len(records) == expected == summaries
    ):
        return True, "manifest and campaign completed"
    return False, "manifest reports incomplete or failed campaign"


def notify(reservation: Reservation | None, subject: str, body: str) -> None:
    mail = Path("/usr/bin/mail")
    if not mail.is_file():
        return
    recipient = os.environ.get("USER") or os.getlogin()
    subprocess.run(
        [str(mail), "-s", subject, recipient], input=body, text=True, check=False
    )


def _describe(reservation: Reservation) -> str:
    return "\n".join(
        [
            f"reservation={reservation.id}",
            f"branch={reservation.branch}",
            f"commit={reservation.commit_sha}",
            f"state={reservation.state}",
            f"logs={_state_dir() / f'reservation-{reservation.id:05d}'}",
            f"control_dir={reservation.control_dir or ''}",
            f"output_root={reservation.output_root or ''}",
            f"error={reservation.error or ''}",
        ]
    )


def run_reservation(store: QueueStore, reservation: Reservation) -> Reservation:
    worktree = ensure_worktree(reservation)
    stdout_log, stderr_log = _reservation_paths(store, reservation)
    runtime = _runtime_python()
    label = f"queue-{reservation.id:05d}-{uuid4().hex[:8]}"
    command = [
        str(runtime),
        str(worktree / "scripts/cs10/parallel_tmux.py"),
        "run",
        "--config",
        str(worktree / reservation.config),
        "--label",
        label,
    ]
    environment = os.environ.copy()
    environment["CS10_RUNTIME_PYTHON"] = str(runtime)
    with (
        stdout_log.open("w", encoding="utf-8") as stdout,
        stderr_log.open("w", encoding="utf-8") as stderr,
    ):
        process = subprocess.Popen(
            command,
            cwd=worktree,
            env=environment,
            stdout=stdout,
            stderr=stderr,
            start_new_session=True,
        )
        store.update(
            reservation.id, state="running", worktree=str(worktree), pid=process.pid
        )
        store.event(reservation.id, "started", " ".join(command))
        returncode = process.wait()
    current = store.get(reservation.id)
    assert current is not None
    if current.state == "cancel_requested":
        store.update(reservation.id, state="cancelled", exit_code=returncode, pid=None)
        current = store.get(reservation.id)
        assert current is not None
        store.event(reservation.id, "cancelled", "process exited after cancellation")
        notify(current, "cs10 queue: job cancelled", _describe(current))
        return current
    try:
        record = _load_launch_record(stdout_log)
        ok, detail = _manifest_succeeded(record)
        store.update(
            reservation.id,
            control_dir=str(record["control_dir"]),
            output_root=str(record["output_root"]),
            exit_code=returncode,
            pid=None,
        )
    except (OSError, ValueError, KeyError, json.JSONDecodeError, RuntimeError) as exc:
        ok, detail = False, str(exc)
    if returncode == 0 and ok:
        store.update(reservation.id, state="succeeded")
        current = store.get(reservation.id)
        assert current is not None
        store.event(reservation.id, "succeeded", detail)
        notify(current, "cs10 queue: job succeeded", _describe(current))
        return current
    store.update(reservation.id, state="failed", error=detail)
    store.set_paused(True)
    current = store.get(reservation.id)
    assert current is not None
    store.event(reservation.id, "failed", detail)
    notify(current, "cs10 queue: job failed; queue paused", _describe(current))
    return current


@contextmanager
def dispatcher_lock(store: QueueStore) -> Iterator[None]:
    lock_path = store.state_dir / "dispatcher.lock"
    with lock_path.open("a+", encoding="utf-8") as lock:
        try:
            fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as exc:
            raise RuntimeError(
                "another cs10 queue dispatcher is already running"
            ) from exc
        yield


def reconcile(store: QueueStore) -> None:
    for reservation in store.list():
        if reservation.state != "running":
            continue
        alive = reservation.pid is not None
        if alive:
            try:
                os.kill(reservation.pid, 0)
            except OSError:
                alive = False
        if not alive:
            store.update(
                reservation.id,
                state="blocked",
                error="dispatcher restart requires operator review",
            )
            store.set_paused(True)
            store.event(
                reservation.id, "blocked", "running process cannot be safely reconciled"
            )


def dispatch(store: QueueStore, *, once: bool, poll_seconds: float) -> int:
    with dispatcher_lock(store):
        reconcile(store)
        while True:
            if store.paused():
                return 0
            reservation = store.next()
            if reservation is None:
                if once:
                    notify(
                        None,
                        "cs10 queue: all jobs complete",
                        "No queued reservations remain.",
                    )
                    return 0
                time.sleep(poll_seconds)
                continue
            run_reservation(store, reservation)
            if once:
                return 0


def cancel(store: QueueStore, reservation_id: int) -> Reservation:
    reservation = store.get(reservation_id)
    if reservation is None:
        raise ValueError(f"unknown reservation: {reservation_id}")
    if reservation.state == "queued":
        store.update(reservation_id, state="cancelled")
        result = store.get(reservation_id)
        assert result is not None
        store.event(reservation_id, "cancelled", "queued reservation cancelled")
        return result
    if reservation.state == "running" and reservation.pid is not None:
        os.killpg(reservation.pid, signal.SIGTERM)
        store.update(reservation_id, state="cancel_requested")
        result = store.get(reservation_id)
        assert result is not None
        store.event(reservation_id, "cancel_requested", "SIGTERM sent to process group")
        return result
    raise ValueError(f"cannot cancel reservation in state: {reservation.state}")


def _print(value: Any) -> None:
    if isinstance(value, Reservation):
        value = asdict(value)
    elif isinstance(value, list):
        value = [asdict(item) for item in value]
    print(json.dumps(value, ensure_ascii=False, indent=2))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    enqueue_parser = subparsers.add_parser("enqueue")
    enqueue_parser.add_argument("--branch", required=True)
    enqueue_parser.add_argument("--config", required=True)
    enqueue_parser.add_argument("--priority", type=int, default=0)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--once", action="store_true")
    run_parser.add_argument("--poll-seconds", type=float, default=30.0)
    subparsers.add_parser("status")
    subparsers.add_parser("pause")
    subparsers.add_parser("resume")
    cancel_parser = subparsers.add_parser("cancel")
    cancel_parser.add_argument("reservation_id", type=int)
    args = parser.parse_args(argv)
    store = QueueStore()
    try:
        if args.action == "enqueue":
            reservation = store.add(
                branch=args.branch,
                commit_sha=resolve_commit(args.branch),
                config=validate_config(args.config),
                priority=args.priority,
            )
            _print(reservation)
        elif args.action == "run":
            return dispatch(store, once=args.once, poll_seconds=args.poll_seconds)
        elif args.action == "status":
            _print(
                {
                    "paused": store.paused(),
                    "reservations": [asdict(item) for item in store.list()],
                }
            )
        elif args.action == "pause":
            store.set_paused(True)
        elif args.action == "resume":
            store.set_paused(False)
        else:
            _print(cancel(store, args.reservation_id))
    finally:
        store.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
