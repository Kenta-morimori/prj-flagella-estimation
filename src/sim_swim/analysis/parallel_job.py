"""Config-level parallel launcher for independent Phase 2 sweep profiles."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import socket
import subprocess
import sys
import time
from typing import Any, Callable, Literal
from uuid import uuid4
from zoneinfo import ZoneInfo

import yaml

from sim_swim.analysis.cli_profiles import load_profile_entry, validate_profile_role


REPO_ROOT = Path(__file__).resolve().parents[3]
SWEEP_DIRECTORY = REPO_ROOT / "conf" / "phase2_sweeps"
THREAD_ENV_KEYS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
WORKER_POLICIES = {"host_cpu", "cs10_qualified"}
SUPPORTED_KINDS = {
    "bundling_alignment",
    "hook_overstretch",
    "motor_scale",
    "shape_stability_grid",
    "single_flagellum_torque",
    "stage_a_2015",
}


@dataclass(frozen=True)
class ParallelJob:
    schema_version: int
    job_id: str
    job_name: str
    config_path: Path
    configs: tuple[Path, ...]
    max_workers: int | Literal["auto"]
    worker_policy: str


@dataclass(frozen=True)
class ResolvedExecution:
    max_workers: int
    worker_policy: str
    thread_environment: dict[str, str]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _now() -> str:
    return datetime.now(ZoneInfo("Asia/Tokyo")).isoformat()


def _relative_to_repo(path: Path) -> str:
    return path.resolve().relative_to(REPO_ROOT).as_posix()


def _require_mapping(value: Any, *, name: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError(f"{name} must be a mapping")
    return value


def _validate_config_path(value: Any) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise ValueError("configs entries must be non-empty paths")
    path = (REPO_ROOT / value).resolve()
    try:
        path.relative_to(SWEEP_DIRECTORY.resolve())
    except ValueError as exc:
        raise ValueError("configs entries must be under conf/phase2_sweeps/") from exc
    if path.suffix not in {".yaml", ".yml"} or not path.is_file():
        raise ValueError(f"sweep config does not exist: {value}")
    entry = load_profile_entry(path)
    validate_profile_role(entry, "sweep")
    if entry["kind"] not in SUPPORTED_KINDS:
        raise ValueError(
            f"unsupported sweep kind {entry['kind']!r}: {_relative_to_repo(path)}"
        )
    return path


def load_parallel_job(path: Path) -> ParallelJob:
    """Load and validate one ``conf/phase2_parallel/<name>/job.yaml`` file."""

    config_path = path.resolve()
    if not config_path.is_file():
        raise ValueError(f"job config does not exist: {path}")
    expected_parent = (REPO_ROOT / "conf" / "phase2_parallel").resolve()
    try:
        job_name = config_path.parent.relative_to(expected_parent).as_posix()
    except ValueError as exc:
        raise ValueError("job config must be under conf/phase2_parallel/") from exc
    if not job_name or "/" in job_name:
        raise ValueError(
            "job config must be directly under conf/phase2_parallel/<job_name>/"
        )

    raw = yaml.safe_load(config_path.read_text(encoding="utf-8")) or {}
    data = _require_mapping(raw, name="job config")
    if data.get("schema_version") != 1:
        raise ValueError("schema_version must be 1")
    job_id = data.get("job_id")
    if not isinstance(job_id, str) or not job_id.strip():
        raise ValueError("job_id must be a non-empty string")
    raw_configs = data.get("configs")
    if not isinstance(raw_configs, list) or not raw_configs:
        raise ValueError("configs must be a non-empty list")
    configs = tuple(_validate_config_path(value) for value in raw_configs)
    if len(set(configs)) != len(configs):
        raise ValueError("configs must not contain duplicates")

    execution = _require_mapping(data.get("execution"), name="execution")
    max_workers = execution.get("max_workers")
    if max_workers != "auto" and (
        not isinstance(max_workers, int)
        or isinstance(max_workers, bool)
        or max_workers < 1
    ):
        raise ValueError("execution.max_workers must be a positive integer or 'auto'")
    worker_policy = execution.get("worker_policy", "host_cpu")
    if worker_policy not in WORKER_POLICIES:
        raise ValueError(
            "execution.worker_policy must be one of: "
            + ", ".join(sorted(WORKER_POLICIES))
        )
    return ParallelJob(
        schema_version=1,
        job_id=job_id.strip(),
        job_name=job_name,
        config_path=config_path,
        configs=configs,
        max_workers=max_workers,
        worker_policy=worker_policy,
    )


def resolve_execution(
    job: ParallelJob, max_workers_override: int | str | None
) -> ResolvedExecution:
    requested = (
        job.max_workers if max_workers_override is None else max_workers_override
    )
    if requested == "auto":
        limit = 8 if job.worker_policy == "cs10_qualified" else (os.cpu_count() or 1)
    elif (
        isinstance(requested, int) and not isinstance(requested, bool) and requested > 0
    ):
        limit = requested
    else:
        raise ValueError("max_workers override must be a positive integer or 'auto'")
    environment = (
        {key: "1" for key in THREAD_ENV_KEYS}
        if job.worker_policy == "cs10_qualified"
        else {}
    )
    return ResolvedExecution(
        max_workers=min(len(job.configs), limit),
        worker_policy=job.worker_policy,
        thread_environment=environment,
    )


def job_output_root(job: ParallelJob, *, output_base_dir: Path | None = None) -> Path:
    base = output_base_dir or REPO_ROOT / "outputs" / "parallel"
    safe_name = job.job_name.replace("/", "_")
    return base / f"{safe_name}__{uuid4().hex[:12]}"


def command_for_config(config: Path, output_dir: Path) -> list[str]:
    entry = load_profile_entry(config)
    output_key = "output_base_dir" if entry["kind"] == "stage_a_2015" else "output_dir"
    return [
        sys.executable,
        "scripts/01_simulate_swimming/run_sweep.py",
        f"config={_relative_to_repo(config)}",
        f"{output_key}={output_dir.resolve()}",
    ]


def build_plan(
    job: ParallelJob, execution: ResolvedExecution, root: Path
) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    for index, config in enumerate(job.configs, start=1):
        config_root = root / f"{index:03d}_{config.stem}"
        run_dir = config_root / "run"
        records.append(
            {
                "index": index,
                "config": _relative_to_repo(config),
                "config_sha256": _sha256(config),
                "kind": load_profile_entry(config)["kind"],
                "command": command_for_config(config, run_dir),
                "output_dir": str(run_dir),
                "stdout_log": str(config_root / "stdout.log"),
                "stderr_log": str(config_root / "stderr.log"),
                "status": "pending",
            }
        )
    return {
        "schema_version": 1,
        "status": "planned",
        "job_id": job.job_id,
        "job_name": job.job_name,
        "job_config": _relative_to_repo(job.config_path),
        "job_config_sha256": _sha256(job.config_path),
        "output_root": str(root),
        "execution": {
            "max_workers": execution.max_workers,
            "worker_policy": execution.worker_policy,
            "thread_environment": {
                key: execution.thread_environment.get(key, os.environ.get(key))
                for key in THREAD_ENV_KEYS
            },
            "thread_environment_overrides": execution.thread_environment,
        },
        "provenance": {
            "hostname": socket.gethostname(),
            "python_executable": sys.executable,
            "python_version": sys.version.split()[0],
        },
        "configs": records,
        "dispatch_order": [],
        "completion_order": [],
        "failed_configs": [],
    }


def _write_manifest(root: Path, manifest: dict[str, Any]) -> None:
    (root / "job_manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )


def _git_provenance() -> dict[str, str | None]:
    result = subprocess.run(
        ["git", "status", "--porcelain=v1", "--branch"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    return {"commit": commit.stdout.strip() or None, "status": result.stdout.strip()}


def run_parallel_job(
    job: ParallelJob,
    execution: ResolvedExecution,
    *,
    output_base_dir: Path | None = None,
    popen: Callable[..., Any] = subprocess.Popen,
    poll_interval_s: float = 0.02,
) -> dict[str, Any]:
    """Run independent profiles, retaining completed results after failures."""

    root = job_output_root(job, output_base_dir=output_base_dir)
    root.mkdir(parents=True, exist_ok=False)
    manifest = build_plan(job, execution, root)
    manifest["status"] = "running"
    manifest["started_at"] = _now()
    manifest["provenance"]["git"] = _git_provenance()
    (root / "job.log").write_text("parallel job started\n", encoding="utf-8")
    _write_manifest(root, manifest)

    environment = os.environ.copy()
    environment.update(execution.thread_environment)
    pending = list(manifest["configs"])
    running: list[tuple[dict[str, Any], Any, Any, Any, float]] = []
    while pending or running:
        while pending and len(running) < execution.max_workers:
            record = pending.pop(0)
            stdout_path = Path(record["stdout_log"])
            stderr_path = Path(record["stderr_log"])
            stdout_path.parent.mkdir(parents=True, exist_ok=False)
            stdout_handle = stdout_path.open("w", encoding="utf-8")
            stderr_handle = stderr_path.open("w", encoding="utf-8")
            record["status"] = "running"
            record["started_at"] = _now()
            started_monotonic = time.monotonic()
            try:
                process = popen(
                    record["command"],
                    cwd=REPO_ROOT,
                    env=environment,
                    stdout=stdout_handle,
                    stderr=stderr_handle,
                )
            except OSError as exc:
                stdout_handle.close()
                stderr_handle.close()
                record["ended_at"] = _now()
                record["wall_time_s"] = time.monotonic() - started_monotonic
                record["exit_code"] = None
                record["status"] = "failed_to_start"
                record["error"] = str(exc)
                manifest["dispatch_order"].append(record["index"])
                manifest["completion_order"].append(record["index"])
                manifest["failed_configs"].append(record["index"])
                _write_manifest(root, manifest)
                continue
            running.append(
                (record, process, stdout_handle, stderr_handle, started_monotonic)
            )
            manifest["dispatch_order"].append(record["index"])
            _write_manifest(root, manifest)
        completed = []
        for item in running:
            record, process, stdout_handle, stderr_handle, started_monotonic = item
            exit_code = process.poll()
            if exit_code is None:
                continue
            stdout_handle.close()
            stderr_handle.close()
            record["ended_at"] = _now()
            record["wall_time_s"] = time.monotonic() - started_monotonic
            record["exit_code"] = exit_code
            record["status"] = "succeeded" if exit_code == 0 else "failed"
            if exit_code != 0:
                manifest["failed_configs"].append(record["index"])
            manifest["completion_order"].append(record["index"])
            completed.append(item)
            _write_manifest(root, manifest)
        running = [item for item in running if item not in completed]
        if running:
            time.sleep(poll_interval_s)
    manifest["ended_at"] = _now()
    manifest["status"] = "succeeded" if not manifest["failed_configs"] else "failed"
    with (root / "job.log").open("a", encoding="utf-8") as handle:
        handle.write(f"parallel job {manifest['status']}\n")
    _write_manifest(root, manifest)
    return manifest
