from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

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
    assert queue.validate_config("conf/job.yaml") == "conf/job.yaml"
    with pytest.raises(ValueError, match="repository-relative"):
        queue.validate_config("../outside.yaml")
    with pytest.raises(ValueError, match="YAML"):
        queue.validate_config("conf/job.txt")


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
                "configs": [{}],
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


def test_cancel_queued_and_running_reservations(monkeypatch, tmp_path: Path) -> None:
    queue = _load_queue("cs10_queue_cancel")
    store = queue.QueueStore(tmp_path)
    try:
        queued = store.add(
            branch="a", commit_sha="a" * 40, config="conf/a.yaml", priority=0
        )
        assert queue.cancel(store, queued.id).state == "cancelled"
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


def test_notify_uses_cs10_mail_for_login_user(monkeypatch, tmp_path: Path) -> None:
    queue = _load_queue("cs10_queue_notify")
    mail = tmp_path / "mail"
    mail.touch()
    monkeypatch.setattr(
        queue, "Path", lambda value: mail if value == "/usr/bin/mail" else Path(value)
    )
    monkeypatch.setenv("USER", "Ktakemori")
    calls: list[tuple[list[str], str]] = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs["input"]))
        return None

    monkeypatch.setattr(queue.subprocess, "run", fake_run)

    queue.notify(None, "subject", "body")

    assert calls == [([str(mail), "-s", "subject", "Ktakemori"], "body")]
