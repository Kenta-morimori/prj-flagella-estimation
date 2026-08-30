#!/usr/bin/env python3
"""Start and inspect cs10-qualified parallel jobs in a persistent tmux session."""

from __future__ import annotations

import argparse
from datetime import datetime
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import tempfile
from typing import Any
from uuid import uuid4
from zoneinfo import ZoneInfo


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = REPOSITORY_ROOT / "src"
CS10_NAS_MOUNT_ROOT = Path("/net/fs01/volume1/work01")
CS10_OUTPUT_BASE = Path(
    "/net/fs01/volume1/work01/Ktakemori/prj-flagella-estimation/outputs"
)
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from sim_swim.analysis.parallel_job import load_parallel_job, resolve_execution  # noqa: E402


def _now() -> datetime:
    return datetime.now(ZoneInfo("Asia/Tokyo"))


def _write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _tmux_path() -> str:
    return shutil.which("tmux") or str(Path.home() / ".local/bin/tmux")


def _require_tmux() -> str:
    path = _tmux_path()
    if not Path(path).is_file():
        raise RuntimeError(
            "tmux is unavailable; install it or add ~/.local/bin to PATH"
        )
    return path


def _git_info() -> dict[str, str]:
    def run(*args: str) -> str:
        return subprocess.run(
            ["git", *args],
            cwd=REPOSITORY_ROOT,
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()

    status = run("status", "--porcelain")
    if status:
        raise RuntimeError("clean worktree is required before cs10 launch")
    return {
        "commit": run("rev-parse", "HEAD"),
        "branch": run("branch", "--show-current"),
    }


def _runtime_python() -> Path:
    configured = os.environ.get("CS10_RUNTIME_PYTHON")
    python = (
        Path(configured) if configured else REPOSITORY_ROOT / ".venv-cs10/bin/python"
    )
    if not python.is_file():
        raise RuntimeError(
            "missing .venv-cs10; run scripts/cs10/setup_environment.sh first"
        )
    environment = os.environ.copy()
    source_path = str(SOURCE_ROOT)
    environment["PYTHONPATH"] = (
        source_path
        if not environment.get("PYTHONPATH")
        else source_path + os.pathsep + environment["PYTHONPATH"]
    )
    subprocess.run(
        [str(python), "-c", "import matplotlib, numpy, yaml, sim_swim"],
        cwd=REPOSITORY_ROOT,
        env=environment,
        check=True,
    )
    return python


def _runtime_environment() -> dict[str, str]:
    """Build the child environment for this worktree's source tree."""
    environment = os.environ.copy()
    source_path = str(SOURCE_ROOT)
    environment["PYTHONPATH"] = (
        source_path
        if not environment.get("PYTHONPATH")
        else source_path + os.pathsep + environment["PYTHONPATH"]
    )
    return environment


def _require_nas_mount() -> None:
    """Reject a path that looks like NAS storage but is not the NAS mount."""
    if not CS10_NAS_MOUNT_ROOT.is_dir() or not CS10_NAS_MOUNT_ROOT.is_mount():
        raise RuntimeError(
            "cs10 NAS mount is unavailable: "
            f"{CS10_NAS_MOUNT_ROOT}; do not create output directories outside "
            "the mounted NAS"
        )


def _require_output_base() -> Path:
    """Return the writable NAS root used for large cs10 job artifacts."""
    _require_nas_mount()
    try:
        CS10_OUTPUT_BASE.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        raise RuntimeError(
            f"could not create cs10 NAS output directory: {CS10_OUTPUT_BASE}"
        ) from exc
    if not CS10_OUTPUT_BASE.is_dir():
        raise RuntimeError(
            "cs10 NAS output directory is unavailable: "
            f"{CS10_OUTPUT_BASE}; ask the cs10 administrators to confirm the mount "
            "and user directory"
        )
    try:
        with tempfile.TemporaryFile(
            mode="w", encoding="utf-8", dir=CS10_OUTPUT_BASE, prefix=".cs10_probe_"
        ) as probe:
            probe.write("ok\n")
            probe.flush()
    except OSError as exc:
        raise RuntimeError(
            f"cs10 NAS output directory is not writable: {CS10_OUTPUT_BASE}"
        ) from exc
    return CS10_OUTPUT_BASE


def _output_paths(label: str, output_base: Path) -> tuple[Path, Path]:
    now = _now()
    control_base = (
        REPOSITORY_ROOT / "outputs" / now.strftime("%Y-%m-%d") / now.strftime("%H%M%S")
    )
    output_timestamp_base = (
        output_base / now.strftime("%Y-%m-%d") / now.strftime("%H%M%S")
    )
    control = control_base / "cs10_parallel" / label
    root = output_timestamp_base / "parallel" / f"{label}__{uuid4().hex[:12]}"
    return control, root


def _prepare_launch(config: Path, *, label: str, session: str | None) -> dict[str, Any]:
    """Validate and record one launch without deciding how it is executed."""
    runtime_python = _runtime_python()
    git = _git_info()
    job = load_parallel_job(config)
    execution = resolve_execution(job, None)
    if execution.worker_policy != "cs10_qualified":
        raise RuntimeError(
            "cs10 tmux helper requires execution.worker_policy=cs10_qualified"
        )
    output_base = _require_output_base()
    control, output_root = _output_paths(label, output_base)
    control.mkdir(parents=True, exist_ok=False)
    command = [
        str(runtime_python),
        "scripts/01_simulate_swimming/run_parallel.py",
        f"config={config.relative_to(REPOSITORY_ROOT)}",
        "--output-root",
        str(output_root),
    ]
    record = {
        "schema_version": 1,
        "status": "started",
        "created_at": _now().isoformat(),
        "session": session,
        "label": label,
        "job_config": str(config.relative_to(REPOSITORY_ROOT)),
        "job_id": job.job_id,
        "output_root": str(output_root),
        "control_dir": str(control),
        "command": command,
        "execution": {
            "max_workers": execution.max_workers,
            "worker_policy": execution.worker_policy,
            "thread_environment": execution.thread_environment,
        },
        "git": git,
    }
    _write_json(control / "launch.json", record)
    return record


def start(config: Path, *, session: str, label: str) -> dict[str, Any]:
    tmux = _require_tmux()
    record = _prepare_launch(config, label=label, session=session)
    existing = subprocess.run([tmux, "has-session", "-t", session], check=False)
    if existing.returncode == 0:
        raise RuntimeError(f"tmux session already exists: {session}")
    control = Path(str(record["control_dir"]))
    output_root = Path(str(record["output_root"]))
    command = list(record["command"])
    launch_script = control / "launch.sh"
    launch_script.write_text(
        "#!/usr/bin/env bash\n"
        "set -u\n"
        f"cd {shlex.quote(str(REPOSITORY_ROOT))}\n"
        + shlex.join(command)
        + f" > {shlex.quote(str(control / 'launcher.stdout.log'))}"
        + f" 2> {shlex.quote(str(control / 'launcher.stderr.log'))}\n"
        "status=$?\n"
        f'printf \'{{"exit_code": %s, "job_root": "%s"}}\\n\' "$status" {shlex.quote(str(output_root))} > {shlex.quote(str(control / "user_exit_marker.json"))}\n'
        f"printf 'Issue parallel launcher finished with exit_code=%s\\n' \"$status\" >> {shlex.quote(str(control / 'launcher.stdout.log'))}\n"
        "exec bash -l\n",
        encoding="utf-8",
    )
    launch_script.chmod(0o700)
    subprocess.run(
        [
            tmux,
            "new-session",
            "-d",
            "-s",
            session,
            "-c",
            str(REPOSITORY_ROOT),
            "bash",
            str(launch_script),
        ],
        check=True,
    )
    return record


def run_foreground(config: Path, *, label: str) -> tuple[dict[str, Any], int]:
    """Run one qualified job in the caller's process tree for queue dispatchers."""
    record = _prepare_launch(config, label=label, session=None)
    control = Path(str(record["control_dir"]))
    with (control / "launcher.stdout.log").open("w", encoding="utf-8") as stdout:
        with (control / "launcher.stderr.log").open("w", encoding="utf-8") as stderr:
            completed = subprocess.run(
                list(record["command"]),
                cwd=REPOSITORY_ROOT,
                env=_runtime_environment(),
                stdout=stdout,
                stderr=stderr,
                check=False,
            )
    _write_json(
        control / "user_exit_marker.json",
        {"exit_code": completed.returncode, "job_root": record["output_root"]},
    )
    return record, completed.returncode


def status(control_dir: Path) -> dict[str, Any]:
    launch = _read_json(control_dir / "launch.json")
    output_root = Path(str(launch["output_root"]))
    result: dict[str, Any] = {
        "launch": launch,
        "marker": None,
        "job": None,
        "campaign": None,
    }
    marker = control_dir / "user_exit_marker.json"
    if marker.is_file():
        result["marker"] = _read_json(marker)
    manifest = output_root / "job_manifest.json"
    if manifest.is_file():
        job = _read_json(manifest)
        records = list(job.get("configs", []) or [])
        result["job"] = {
            "status": job.get("status"),
            "failed_configs": job.get("failed_configs", []),
            "dispatched": len(job.get("dispatch_order", [])),
            "completed": len(job.get("completion_order", [])),
            "succeeded": sum(item.get("status") == "succeeded" for item in records),
            "failed": sum(
                str(item.get("status", "")).startswith("failed") for item in records
            ),
            "aggregation": job.get("aggregation"),
        }
    completion = output_root / "campaign/campaign_completion.json"
    conditions = output_root / "campaign/conditions"
    if completion.is_file():
        result["campaign"] = {
            **_read_json(completion),
            "run_summary_count": (
                sum(
                    (condition / "run_summary.json").is_file()
                    for condition in conditions.iterdir()
                )
                if conditions.is_dir()
                else 0
            ),
        }
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="action", required=True)
    start_parser = subparsers.add_parser("start")
    start_parser.add_argument("--config", type=Path, required=True)
    start_parser.add_argument("--session", default="issue203")
    start_parser.add_argument("--label", required=True)
    run_parser = subparsers.add_parser("run")
    run_parser.add_argument("--config", type=Path, required=True)
    run_parser.add_argument("--label", required=True)
    status_parser = subparsers.add_parser("status")
    status_parser.add_argument("--control-dir", type=Path, required=True)
    attach_parser = subparsers.add_parser("attach")
    attach_parser.add_argument("--session", required=True)
    args = parser.parse_args(argv)
    if args.action == "start":
        print(
            json.dumps(
                start(args.config.resolve(), session=args.session, label=args.label),
                ensure_ascii=False,
                indent=2,
            )
        )
        return 0
    if args.action == "status":
        print(
            json.dumps(status(args.control_dir.resolve()), ensure_ascii=False, indent=2)
        )
        return 0
    if args.action == "run":
        record, returncode = run_foreground(args.config.resolve(), label=args.label)
        print(json.dumps(record, ensure_ascii=False, indent=2))
        return returncode
    tmux = _require_tmux()
    return subprocess.run([tmux, "attach", "-t", args.session], check=False).returncode


if __name__ == "__main__":
    raise SystemExit(main())
