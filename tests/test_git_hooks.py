from __future__ import annotations

import os
import stat
import subprocess
from pathlib import Path


HOOK_PATH = Path(__file__).resolve().parents[1] / ".githooks" / "pre-commit"


def _write_executable(path: Path, content: str) -> None:
    path.write_text(content, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR)


def _make_fake_uv(tmp_path: Path, pytest_exit_code: int) -> Path:
    uv_path = tmp_path / "uv"
    _write_executable(
        uv_path,
        f"""#!/bin/sh
if [ "$1" != "run" ]; then
  exit 99
fi
shift
if [ "$1" = "ruff" ]; then
  exit 0
fi
if [ "$1" = "pytest" ]; then
  exit {pytest_exit_code}
fi
exit 98
""",
    )
    return uv_path


def _make_repo_with_tests(tmp_path: Path) -> Path:
    repo = tmp_path / "repo"
    (repo / "tests").mkdir(parents=True)
    (repo / "tests" / "test_dummy.py").write_text("def test_dummy():\n    pass\n")
    return repo


def _run_pre_commit(repo: Path, uv_dir: Path) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env["PATH"] = f"{uv_dir}{os.pathsep}{env['PATH']}"
    return subprocess.run(
        ["sh", str(HOOK_PATH)],
        cwd=repo,
        env=env,
        text=True,
        capture_output=True,
    )


def test_pre_commit_fails_when_light_marker_collects_no_tests(tmp_path: Path) -> None:
    repo = _make_repo_with_tests(tmp_path)
    uv_path = _make_fake_uv(tmp_path, pytest_exit_code=5)

    result = _run_pre_commit(repo, uv_path.parent)

    assert result.returncode == 5
    assert "no light-marked tests were collected" in result.stderr


def test_pre_commit_succeeds_when_light_marker_collects_tests(tmp_path: Path) -> None:
    repo = _make_repo_with_tests(tmp_path)
    uv_path = _make_fake_uv(tmp_path, pytest_exit_code=0)

    result = _run_pre_commit(repo, uv_path.parent)

    assert result.returncode == 0
