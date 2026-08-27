from __future__ import annotations

import importlib.util
from pathlib import Path
import subprocess
import sys


def _module():
    path = Path(__file__).resolve().parents[1] / "scripts/cs10/queue.py"
    spec = importlib.util.spec_from_file_location("cs10_queue", path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules["cs10_queue"] = module
    spec.loader.exec_module(module)
    return module


def _repo(path: Path) -> None:
    subprocess.run(["git", "init"], cwd=path, check=True, capture_output=True)
    subprocess.run(
        ["git", "config", "user.email", "test@example.com"], cwd=path, check=True
    )
    subprocess.run(["git", "config", "user.name", "Test"], cwd=path, check=True)
    (path / "a").write_text("a", encoding="utf-8")
    subprocess.run(["git", "add", "a"], cwd=path, check=True)
    subprocess.run(
        ["git", "commit", "-m", "initial"], cwd=path, check=True, capture_output=True
    )


def test_enqueue_resolves_branch_to_immutable_commit(tmp_path: Path) -> None:
    queue = _module()
    repo = tmp_path / "repo"
    repo.mkdir()
    _repo(repo)
    db = queue._connect(tmp_path / "state")
    job_id = queue.enqueue(db, repo, "HEAD", "true")
    job = db.execute("SELECT * FROM jobs WHERE id=?", (job_id,)).fetchone()
    assert job["status"] == "queued"
    assert len(job["commit_sha"]) == 40


def test_manifest_success_contract(tmp_path: Path) -> None:
    queue = _module()
    (tmp_path / "job_manifest.json").write_text(
        '{"status":"succeeded","failed_configs":[]}', encoding="utf-8"
    )
    assert queue._manifest_succeeded(tmp_path) == (True, "")
    (tmp_path / "job_manifest.json").write_text(
        '{"status":"failed","failed_configs":[]}', encoding="utf-8"
    )
    assert not queue._manifest_succeeded(tmp_path)[0]
