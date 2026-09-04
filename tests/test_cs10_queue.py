from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]


def _load_queue(name: str):
    path = ROOT / "scripts/cs10/queue.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def test_store_orders_by_priority_then_fifo_and_persists_pause(tmp_path: Path) -> None:
    queue = _load_queue("cs10_queue_store")
    store = queue.QueueStore(tmp_path)
    try:
        first = store.add(
            branch="origin/a", commit_sha="a" * 40, config="conf/a.yaml", priority=0
        )
        highest = store.add(
            branch="origin/b", commit_sha="b" * 40, config="conf/b.yaml", priority=2
        )
        second = store.add(
            branch="origin/c", commit_sha="c" * 40, config="conf/c.yaml", priority=0
        )
        assert store.next().id == highest.id
        assert [item.id for item in store.list()] == [highest.id, first.id, second.id]
        store.set_paused(True)
    finally:
        store.close()
    reopened = queue.QueueStore(tmp_path)
    try:
        assert reopened.paused()
    finally:
        reopened.close()


def test_validate_config_requires_checked_in_relative_yaml(
    monkeypatch, tmp_path: Path
) -> None:
    queue = _load_queue("cs10_queue_config")
    repository = tmp_path / "repository"
    config = repository / "conf/job.yaml"
    config.parent.mkdir(parents=True)
    config.write_text("job_id: example\n", encoding="utf-8")
    monkeypatch.setattr(queue, "REPOSITORY_ROOT", repository)
    commands: list[tuple[str, ...]] = []
    monkeypatch.setattr(queue, "_git", lambda *args: commands.append(args) or "")
    assert queue.validate_config("conf/job.yaml", "a" * 40) == "conf/job.yaml"
    assert commands == [("cat-file", "-e", f"{'a' * 40}:conf/job.yaml")]
    with pytest.raises(ValueError, match="repository-relative"):
        queue.validate_config("../outside.yaml", "a" * 40)
    with pytest.raises(ValueError, match="YAML"):
        queue.validate_config("conf/job.txt", "a" * 40)


def test_ensure_worktree_uses_fixed_commit(monkeypatch, tmp_path: Path) -> None:
    queue = _load_queue("cs10_queue_worktree")
    monkeypatch.setattr(queue, "_worktree_dir", lambda: tmp_path / "worktrees")
    commands: list[tuple[str, ...]] = []
    monkeypatch.setattr(queue, "_git", lambda *args: commands.append(args) or "")
    reservation = queue.Reservation(
        7,
        "origin/topic",
        "abc123def4567890",
        "conf/job.yaml",
        0,
        "queued",
        "now",
        None,
        None,
        None,
        None,
        None,
        None,
    )
    path = queue.ensure_worktree(reservation)
    assert path.name == "queue-00007-abc123def456"
    assert commands == [
        ("worktree", "add", "--detach", str(path), reservation.commit_sha)
    ]


def test_manifest_success_requires_all_condition_summaries(tmp_path: Path) -> None:
    queue = _load_queue("cs10_queue_manifest")
    root = tmp_path / "output"
    (root / "campaign/conditions/a").mkdir(parents=True)
    (root / "campaign/conditions/a/run_summary.json").write_text("{}", encoding="utf-8")
    (root / "job_manifest.json").write_text(
        json.dumps(
            {
                "status": "succeeded",
                "failed_configs": [],
                "configs": [{"status": "succeeded"}],
                "aggregation": {"status": "completed"},
            }
        ),
        encoding="utf-8",
    )
    (root / "campaign/campaign_completion.json").write_text(
        json.dumps(
            {"status": "completed", "exit_code": 0, "expected_condition_count": 1}
        ),
        encoding="utf-8",
    )
    assert queue._manifest_succeeded({"output_root": str(root)}) == (
        True,
        "manifest and campaign completed",
    )
    (root / "campaign/conditions/a/run_summary.json").unlink()
    assert queue._manifest_succeeded({"output_root": str(root)})[0] is False


def test_manifest_success_allows_completed_nongeneric_parallel_job(
    tmp_path: Path,
) -> None:
    queue = _load_queue("cs10_queue_nongeneric_manifest")
    root = tmp_path / "output"
    root.mkdir()
    (root / "job_manifest.json").write_text(
        json.dumps(
            {
                "status": "succeeded",
                "failed_configs": [],
                "configs": [{"status": "succeeded"}],
            }
        ),
        encoding="utf-8",
    )
    assert queue._manifest_succeeded({"output_root": str(root)}) == (
        True,
        "job manifest completed",
    )


def test_cancel_queued_and_running_reservations(monkeypatch, tmp_path: Path) -> None:
    queue = _load_queue("cs10_queue_cancel")
    store = queue.QueueStore(tmp_path)
    try:
        notifications: list[int | None] = []
        monkeypatch.setattr(
            queue,
            "_notify",
            lambda _store, reservation, *_: notifications.append(reservation.id),
        )
        queued = store.add(
            branch="a", commit_sha="a" * 40, config="conf/a.yaml", priority=0
        )
        assert queue.cancel(store, queued.id).state == "cancelled"
        assert notifications == [queued.id]
        running = store.add(
            branch="b", commit_sha="b" * 40, config="conf/b.yaml", priority=0
        )
        store.update(running.id, state="running", pid=1234)
        sent: list[tuple[int, int]] = []
        monkeypatch.setattr(
            queue.os, "killpg", lambda pid, sig: sent.append((pid, sig))
        )
        assert queue.cancel(store, running.id).state == "cancel_requested"
        assert sent == [(1234, queue.signal.SIGTERM)]
    finally:
        store.close()


def test_reconcile_blocks_unknown_running_process_and_pauses_queue(
    monkeypatch, tmp_path: Path
) -> None:
    queue = _load_queue("cs10_queue_reconcile")
    store = queue.QueueStore(tmp_path)
    try:
        reservation = store.add(
            branch="a", commit_sha="a" * 40, config="conf/a.yaml", priority=0
        )
        store.update(reservation.id, state="running", pid=9999)
        monkeypatch.setattr(
            queue.os, "kill", lambda *_: (_ for _ in ()).throw(OSError())
        )
        queue.reconcile(store)
        assert store.get(reservation.id).state == "blocked"
        assert store.paused()
    finally:
        store.close()


def test_reconcile_blocks_live_running_process_and_pauses_queue(
    monkeypatch, tmp_path: Path
) -> None:
    queue = _load_queue("cs10_queue_live_reconcile")
    store = queue.QueueStore(tmp_path)
    try:
        reservation = store.add(
            branch="a", commit_sha="a" * 40, config="conf/a.yaml", priority=0
        )
        store.update(reservation.id, state="running", pid=9999)
        monkeypatch.setattr(queue.os, "kill", lambda *_: None)
        queue.reconcile(store)
        result = store.get(reservation.id)
        assert result is not None
        assert result.state == "blocked"
        assert "live child" in str(result.error)
        assert store.paused()
    finally:
        store.close()


def test_dispatcher_lock_rejects_a_second_dispatcher(tmp_path: Path) -> None:
    queue = _load_queue("cs10_queue_lock")
    store = queue.QueueStore(tmp_path)
    try:
        with queue.dispatcher_lock(store):
            with pytest.raises(RuntimeError, match="already running"):
                with queue.dispatcher_lock(store):
                    pass
    finally:
        store.close()


def test_notification_setup_requires_github_cli_and_authentication(monkeypatch) -> None:
    queue = _load_queue("cs10_queue_notify_setup")
    monkeypatch.setattr(queue.shutil, "which", lambda _: None)
    with pytest.raises(RuntimeError, match="GitHub CLI is unavailable"):
        queue.validate_notification_setup()

    monkeypatch.setattr(queue.shutil, "which", lambda _: "/usr/bin/gh")
    monkeypatch.setattr(
        queue.subprocess,
        "run",
        lambda command, **kwargs: subprocess.CompletedProcess(command, 1, "", ""),
    )
    with pytest.raises(RuntimeError, match="not authenticated"):
        queue.validate_notification_setup()


def test_notify_dispatches_actions_workflow_without_recipient(
    monkeypatch, tmp_path: Path
) -> None:
    queue = _load_queue("cs10_queue_notify_dispatch")
    monkeypatch.setattr(queue, "validate_notification_setup", lambda: None)
    monkeypatch.setattr(queue, "_state_dir", lambda: tmp_path / "state")
    commands: list[list[str]] = []

    def fake_run(command, **kwargs):
        commands.append(command)
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(queue.subprocess, "run", fake_run)
    reservation = queue.Reservation(
        7,
        "origin/topic",
        "a" * 40,
        "conf/job.yaml",
        0,
        "failed",
        "now",
        None,
        "/control",
        "/output",
        None,
        1,
        "manifest failed",
    )
    queue.notify(reservation, "failed")

    assert commands == [
        [
            "gh",
            "workflow",
            "run",
            "cs10-queue-notify.yml",
            "--repo",
            "Kenta-morimori/prj-flagella-estimation",
            "--ref",
            "main",
            "--field",
            "event_type=failed",
            "--field",
            "reservation_id=7",
            "--field",
            "branch=origin/topic",
            "--field",
            f"commit_sha={'a' * 40}",
            "--field",
            "state=failed",
            "--field",
            f"logs={tmp_path}/state/reservation-00007",
            "--field",
            "control_dir=/control",
            "--field",
            "output_root=/output",
            "--field",
            "error=manifest failed",
        ]
    ]
    assert not any("@" in value for value in commands[0])

    monkeypatch.setattr(
        queue.subprocess,
        "run",
        lambda command, **kwargs: subprocess.CompletedProcess(command, 1, "", ""),
    )
    with pytest.raises(RuntimeError, match="dispatch failed with exit code 1"):
        queue.notify(reservation, "failed")


def test_notify_records_dispatch_success_and_failure(
    monkeypatch, tmp_path: Path
) -> None:
    queue = _load_queue("cs10_queue_notify_events")
    store = queue.QueueStore(tmp_path)
    try:
        monkeypatch.setattr(queue, "notify", lambda *args: None)
        queue._notify(store, None, "queue_completed")
        monkeypatch.setattr(
            queue,
            "notify",
            lambda *args: (_ for _ in ()).throw(RuntimeError("dispatch failed")),
        )
        queue._notify(store, None, "queue_completed")
        events = store.connection.execute(
            "SELECT kind, detail FROM events ORDER BY id"
        ).fetchall()
        assert [(row["kind"], row["detail"]) for row in events] == [
            (
                "notification_dispatched",
                "GitHub Actions accepted queue_completed notification dispatch",
            ),
            ("notification_failed", "dispatch failed"),
        ]
    finally:
        store.close()


def test_notify_records_oserror_without_recipient(monkeypatch, tmp_path: Path) -> None:
    queue = _load_queue("cs10_queue_notify_config_oserror")
    store = queue.QueueStore(tmp_path)
    try:
        monkeypatch.setattr(
            queue,
            "validate_notification_setup",
            lambda: (_ for _ in ()).throw(OSError("authentication cache unreadable")),
        )
        queue._notify(store, None, "queue_completed")
        event = store.connection.execute(
            "SELECT kind, detail FROM events ORDER BY id"
        ).fetchone()
        assert event["kind"] == "notification_failed"
        assert "OSError" in event["detail"]
    finally:
        store.close()


def test_dispatch_runs_two_reservations_sequentially(
    monkeypatch, tmp_path: Path
) -> None:
    queue = _load_queue("cs10_queue_sequential_dispatch")
    store = queue.QueueStore(tmp_path)
    try:
        first = store.add(
            branch="a", commit_sha="a" * 40, config="conf/a.yaml", priority=0
        )
        second = store.add(
            branch="b", commit_sha="b" * 40, config="conf/b.yaml", priority=0
        )
        started: list[int] = []
        monkeypatch.setattr(queue, "validate_notification_setup", lambda: None)
        monkeypatch.setattr(queue, "_notify", lambda *args: None)

        def fake_run_reservation(current_store, reservation):
            started.append(reservation.id)
            current_store.update(reservation.id, state="succeeded")
            result = current_store.get(reservation.id)
            assert result is not None
            return result

        monkeypatch.setattr(queue, "run_reservation", fake_run_reservation)

        class DispatcherFinished(Exception):
            pass

        monkeypatch.setattr(
            queue.time,
            "sleep",
            lambda _: (_ for _ in ()).throw(DispatcherFinished()),
        )
        with pytest.raises(DispatcherFinished):
            queue.dispatch(store, once=False, poll_seconds=0)

        assert started == [first.id, second.id]
        assert [item.state for item in store.list()] == ["succeeded", "succeeded"]
    finally:
        store.close()
