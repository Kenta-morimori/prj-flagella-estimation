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
import shutil
import sqlite3
import subprocess
import time
from typing import Any, Iterator
from uuid import uuid4
from zoneinfo import ZoneInfo


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_STATE_DIR = Path.home() / ".local/state/prj-flagella-estimation/cs10-queue"
DEFAULT_WORKTREE_DIR = Path.home() / "src/prj-flagella-estimation-queue-worktrees"
GITHUB_BINARY = "gh"
NOTIFICATION_REPOSITORY = "Kenta-morimori/prj-flagella-estimation"
NOTIFICATION_WORKFLOW = "cs10-queue-notify.yml"
NOTIFICATION_REF = "main"
TERMINAL_STATES = {"succeeded", "failed", "cancelled", "blocked"}
NOTIFIABLE_STATES = {"succeeded", "failed", "cancelled"}


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
    notification_attempted: bool = False


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
                error TEXT,
                notification_attempted INTEGER NOT NULL DEFAULT 0
            )"""
        )
        columns = {
            row["name"]
            for row in self.connection.execute("PRAGMA table_info(reservations)")
        }
        if "notification_attempted" not in columns:
            self.connection.execute(
                "ALTER TABLE reservations ADD COLUMN notification_attempted INTEGER NOT NULL DEFAULT 0"
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

    def claim_notification(self, reservation_id: int) -> bool:
        """Persist a single automatic notification attempt for one reservation."""
        cursor = self.connection.execute(
            "UPDATE reservations SET notification_attempted=1 "
            "WHERE id=? AND notification_attempted=0",
            (reservation_id,),
        )
        self.connection.commit()
        return cursor.rowcount == 1

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


def validate_config(config: str, commit_sha: str) -> str:
    path = Path(config)
    if path.is_absolute() or ".." in path.parts or path.suffix not in {".yaml", ".yml"}:
        raise ValueError("config must be a repository-relative YAML path")
    if REPOSITORY_ROOT not in (REPOSITORY_ROOT / path).resolve().parents:
        raise ValueError(f"parallel-job config is unavailable: {config}")
    try:
        _git("cat-file", "-e", f"{commit_sha}:{path}")
    except subprocess.CalledProcessError as exc:
        raise ValueError(
            f"parallel-job config is unavailable at {commit_sha}: {config}"
        ) from exc
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
    if not manifest.is_file():
        return False, "missing job manifest"
    job = json.loads(manifest.read_text(encoding="utf-8"))
    failed = job.get("failed_configs", [])
    records = list(job.get("configs", []) or [])
    if (
        job.get("status") != "succeeded"
        or failed
        or not records
        or any(record.get("status") != "succeeded" for record in records)
    ):
        return False, "job manifest reports incomplete or failed configs"
    if not completion.is_file():
        return True, "job manifest completed"
    campaign = json.loads(completion.read_text(encoding="utf-8"))
    expected = int(campaign.get("expected_condition_count", -1))
    conditions = root / "campaign/conditions"
    summaries = (
        sum((item / "run_summary.json").is_file() for item in conditions.iterdir())
        if conditions.is_dir()
        else 0
    )
    if (
        job.get("aggregation", {}).get("status") == "completed"
        and campaign.get("status") == "completed"
        and campaign.get("exit_code") == 0
        and len(records) == expected == summaries
    ):
        return True, "manifest and campaign completed"
    return False, "manifest reports incomplete or failed campaign"


def validate_notification_setup() -> None:
    """Reject dispatcher startup unless it can dispatch the Actions workflow."""
    binary = shutil.which(GITHUB_BINARY)
    if binary is None:
        raise RuntimeError(f"required GitHub CLI is unavailable: {GITHUB_BINARY}")
    result = subprocess.run(
        [binary, "auth", "status", "--hostname", "github.com"],
        text=True,
        capture_output=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError("GitHub CLI is not authenticated for github.com")


def notify(reservation: Reservation) -> None:
    """Dispatch one final reservation notification without exposing its recipient."""
    try:
        validate_notification_setup()
        reservation_log_dir = _state_dir() / f"reservation-{reservation.id:05d}"
        result = subprocess.run(
            [
                GITHUB_BINARY,
                "workflow",
                "run",
                NOTIFICATION_WORKFLOW,
                "--repo",
                NOTIFICATION_REPOSITORY,
                "--ref",
                NOTIFICATION_REF,
                "--raw-field",
                "event_type=job_completed",
                "--raw-field",
                f"reservation_id={reservation.id}",
                "--raw-field",
                f"branch={reservation.branch}",
                "--raw-field",
                f"commit_sha={reservation.commit_sha}",
                "--raw-field",
                f"state={reservation.state}",
                "--raw-field",
                f"logs={reservation_log_dir}",
                "--raw-field",
                f"control_dir={reservation.control_dir or ''}",
                "--raw-field",
                f"output_root={reservation.output_root or ''}",
                "--raw-field",
                f"error={reservation.error or ''}",
            ],
            text=True,
            capture_output=True,
            check=False,
        )
    except OSError as exc:
        raise RuntimeError(
            f"GitHub Actions notification could not start: {exc.__class__.__name__}"
        ) from exc
    if result.returncode != 0:
        raise RuntimeError(
            f"GitHub Actions dispatch failed with exit code {result.returncode}"
        )


def _notify_completion(store: QueueStore, reservation: Reservation) -> None:
    """Dispatch at most one final notification after a reservation is terminal."""
    if reservation.state not in NOTIFIABLE_STATES:
        return
    if not store.claim_notification(reservation.id):
        return
    try:
        notify(reservation)
    except RuntimeError as exc:
        store.event(reservation.id, "notification_failed", str(exc))
    else:
        store.event(
            reservation.id,
            "notification_dispatched",
            "GitHub Actions accepted job_completed notification dispatch",
        )


def _finalize_reservation(
    store: QueueStore,
    reservation_id: int,
    *,
    state: str,
    detail: str,
    pause_queue: bool = False,
    **values: Any,
) -> Reservation:
    """Persist a final state, then send its one permitted notification."""
    if state not in NOTIFIABLE_STATES:
        raise ValueError(f"state is not a notifiable terminal state: {state}")
    store.update(reservation_id, state=state, **values)
    if pause_queue:
        store.set_paused(True)
    reservation = store.get(reservation_id)
    assert reservation is not None
    store.event(reservation.id, state, detail)
    _notify_completion(store, reservation)
    return reservation


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
        return _finalize_reservation(
            store,
            reservation.id,
            state="cancelled",
            detail="process exited after cancellation",
            exit_code=returncode,
            pid=None,
        )
    final_values: dict[str, Any] = {"exit_code": returncode, "pid": None}
    try:
        record = _load_launch_record(stdout_log)
        ok, detail = _manifest_succeeded(record)
        final_values.update(
            control_dir=str(record["control_dir"]),
            output_root=str(record["output_root"]),
        )
    except (OSError, ValueError, KeyError, json.JSONDecodeError, RuntimeError) as exc:
        ok, detail = False, str(exc)
    if returncode == 0 and ok:
        return _finalize_reservation(
            store, reservation.id, state="succeeded", detail=detail, **final_values
        )
    return _finalize_reservation(
        store,
        reservation.id,
        state="failed",
        detail=detail,
        error=detail,
        pause_queue=True,
        **final_values,
    )


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
        alive = False
        if reservation.pid is not None:
            try:
                os.kill(reservation.pid, 0)
                alive = True
            except OSError:
                pass
        detail = (
            "dispatcher restart found a live child process; operator review required"
            if alive
            else "dispatcher restart found an exited child process; operator review required"
        )
        store.update(reservation.id, state="blocked", error=detail)
        store.set_paused(True)
        store.event(reservation.id, "blocked", detail)


def dispatch(store: QueueStore, *, once: bool, poll_seconds: float) -> int:
    validate_notification_setup()
    with dispatcher_lock(store):
        reconcile(store)
        while True:
            if store.paused():
                return 0
            reservation = store.next()
            if reservation is None:
                if once:
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
        return _finalize_reservation(
            store,
            reservation_id,
            state="cancelled",
            detail="queued reservation cancelled",
        )
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
            commit_sha = resolve_commit(args.branch)
            reservation = store.add(
                branch=args.branch,
                commit_sha=commit_sha,
                config=validate_config(args.config, commit_sha),
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
