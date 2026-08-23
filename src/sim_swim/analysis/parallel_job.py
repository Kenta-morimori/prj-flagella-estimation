"""Config-level parallel launcher for independent Phase 2 sweep profiles."""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
import csv
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
from sim_swim.analysis.multi_run_campaign import (
    apply_campaign_cli_overrides,
    build_campaign_conditions,
    campaign_axes_metadata,
    load_yaml,
)
from sim_swim.analysis.sweeps.generic_multi_run import (
    _condition_row,
    _manifest_condition_record,
    _summary_fieldnames,
)
from sim_swim.sim.params import SimulationConfig


REPO_ROOT = Path(__file__).resolve().parents[3]
SWEEP_DIRECTORY = REPO_ROOT / "conf" / "phase2_sweeps"
MULTI_RUN_DIRECTORY = REPO_ROOT / "conf" / "phase2_multi_run"
THREAD_ENV_KEYS = ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
WORKER_POLICIES = {"host_cpu", "cs10_qualified"}
JOB_KEYS = {"schema_version", "job_id", "configs", "conditions", "execution"}
EXECUTION_KEYS = {"max_workers", "worker_policy"}
CONFIG_ENTRY_KEYS = {"path", "overrides"}
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
    config_overrides: dict[Path, tuple[str, ...]] = field(default_factory=dict)
    condition_ids: tuple[str, ...] = ()

    @property
    def task_count(self) -> int:
        return len(self.condition_ids) if self.condition_ids else len(self.configs)

    @property
    def is_generic_campaign_job(self) -> bool:
        return bool(self.condition_ids)


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


def _reject_unknown_keys(data: dict[str, Any], *, allowed: set[str], name: str) -> None:
    unknown = sorted(set(data) - allowed)
    if unknown:
        raise ValueError(f"{name} contains unsupported keys: {', '.join(unknown)}")


def _validate_config_path(value: Any) -> Path:
    if not isinstance(value, str) or not value.strip():
        raise ValueError("configs entries must be non-empty paths")
    path = (REPO_ROOT / value).resolve()
    is_sweep = path.is_relative_to(SWEEP_DIRECTORY.resolve())
    is_multi_run = path.is_relative_to(MULTI_RUN_DIRECTORY.resolve())
    if not (is_sweep or is_multi_run):
        raise ValueError(
            "configs entries must be under conf/phase2_sweeps/ or conf/phase2_multi_run/"
        )
    if path.suffix not in {".yaml", ".yml"} or not path.is_file():
        raise ValueError(f"config does not exist: {value}")
    entry = load_profile_entry(path)
    validate_profile_role(entry, "sweep")
    if entry["kind"] not in {*SUPPORTED_KINDS, "generic_multi_run"}:
        raise ValueError(
            f"unsupported config kind {entry['kind']!r}: {_relative_to_repo(path)}"
        )
    return path


def _validate_config_overrides(value: Any) -> tuple[str, ...]:
    if value is None:
        return ()
    if not isinstance(value, list):
        raise ValueError("config overrides must be a list")
    overrides: list[str] = []
    keys: set[str] = set()
    for item in value:
        if not isinstance(item, str) or not item.strip():
            raise ValueError("config overrides must contain non-empty strings")
        key, separator, raw_value = item.partition("=")
        if not separator or not key or not raw_value or key.strip() != key:
            raise ValueError("config overrides must use KEY=VALUE")
        normalized_key = key.replace("-", "_")
        if normalized_key in {"output_dir", "output_base_dir"}:
            raise ValueError(
                "config overrides must not set launcher-managed output paths"
            )
        if normalized_key in keys:
            raise ValueError(f"config overrides must not repeat key: {key}")
        keys.add(normalized_key)
        overrides.append(item)
    return tuple(overrides)


def _validate_config_entry(value: Any) -> tuple[Path, tuple[str, ...]]:
    if isinstance(value, str):
        return _validate_config_path(value), ()
    data = _require_mapping(value, name="config entry")
    _reject_unknown_keys(data, allowed=CONFIG_ENTRY_KEYS, name="config entry")
    return _validate_config_path(data.get("path")), _validate_config_overrides(
        data.get("overrides")
    )


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
    _reject_unknown_keys(data, allowed=JOB_KEYS, name="job config")
    if data.get("schema_version") != 1:
        raise ValueError("schema_version must be 1")
    job_id = data.get("job_id")
    if not isinstance(job_id, str) or not job_id.strip():
        raise ValueError("job_id must be a non-empty string")
    raw_configs = data.get("configs")
    if not isinstance(raw_configs, list) or not raw_configs:
        raise ValueError("configs must be a non-empty list")
    config_entries = tuple(_validate_config_entry(value) for value in raw_configs)
    configs = tuple(entry[0] for entry in config_entries)
    if len(set(configs)) != len(configs):
        raise ValueError("configs must not contain duplicates")
    generic_configs = [
        config
        for config in configs
        if load_profile_entry(config)["kind"] == "generic_multi_run"
    ]
    raw_conditions = data.get("conditions")
    if generic_configs:
        if len(configs) != 1:
            raise ValueError(
                "generic_multi_run parallel jobs must contain exactly one config"
            )
        if not isinstance(raw_conditions, list) or not raw_conditions:
            raise ValueError(
                "generic_multi_run parallel jobs require non-empty conditions"
            )
        condition_ids = tuple(str(item) for item in raw_conditions)
        if any(not item.strip() for item in condition_ids):
            raise ValueError("conditions must contain non-empty condition IDs")
        if len(set(condition_ids)) != len(condition_ids):
            raise ValueError("conditions must not contain duplicates")
        effective = apply_campaign_cli_overrides(
            load_yaml(generic_configs[0]), list(config_entries[0][1])
        )
        available = {
            item["condition_id"] for item in build_campaign_conditions(effective)
        }
        unknown = [item for item in condition_ids if item not in available]
        if unknown:
            raise ValueError(
                "conditions contains unknown condition IDs: " + ", ".join(unknown)
            )
    else:
        if raw_conditions is not None:
            raise ValueError("conditions is supported only for generic_multi_run jobs")
        condition_ids = ()

    execution = _require_mapping(data.get("execution"), name="execution")
    _reject_unknown_keys(execution, allowed=EXECUTION_KEYS, name="execution")
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
        config_overrides={config: overrides for config, overrides in config_entries},
        condition_ids=condition_ids,
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
    if job.worker_policy == "cs10_qualified":
        limit = min(limit, 8)
    environment = (
        {key: "1" for key in THREAD_ENV_KEYS}
        if job.worker_policy == "cs10_qualified"
        else {}
    )
    return ResolvedExecution(
        max_workers=min(job.task_count, limit),
        worker_policy=job.worker_policy,
        thread_environment=environment,
    )


def job_output_root(job: ParallelJob, *, output_base_dir: Path | None = None) -> Path:
    base = output_base_dir or REPO_ROOT / "outputs"
    safe_name = job.job_name.replace("/", "_")
    now = datetime.now(ZoneInfo("Asia/Tokyo"))
    return (
        base
        / now.strftime("%Y-%m-%d")
        / now.strftime("%H%M%S")
        / "parallel"
        / f"{safe_name}__{uuid4().hex[:12]}"
    )


def command_for_config(
    config: Path, output_dir: Path, overrides: tuple[str, ...] = ()
) -> list[str]:
    entry = load_profile_entry(config)
    output_key = "output_base_dir" if entry["kind"] == "stage_a_2015" else "output_dir"
    return [
        sys.executable,
        "scripts/01_simulate_swimming/run_sweep.py",
        f"config={_relative_to_repo(config)}",
        f"{output_key}={output_dir.resolve()}",
        *overrides,
    ]


def _generic_command(
    config: Path, output_dir: Path, condition_id: str, overrides: tuple[str, ...]
) -> list[str]:
    return [
        sys.executable,
        "scripts/01_simulate_swimming/run_multi_run.py",
        f"config={_relative_to_repo(config)}",
        f"output.base_dir={output_dir.resolve()}",
        "output.timestamp_subdir=false",
        f"sweep.include_condition_ids=[{condition_id}]",
        *overrides,
    ]


def build_plan(
    job: ParallelJob, execution: ResolvedExecution, root: Path
) -> dict[str, Any]:
    records: list[dict[str, Any]] = []
    task_items = (
        [(job.configs[0], condition_id) for condition_id in job.condition_ids]
        if job.is_generic_campaign_job
        else [(config, None) for config in job.configs]
    )
    for index, (config, condition_id) in enumerate(task_items, start=1):
        suffix = condition_id if condition_id is not None else config.stem
        config_root = root / "children" / f"{index:03d}_{suffix}"
        run_dir = config_root / "run"
        kind = load_profile_entry(config)["kind"]
        overrides = job.config_overrides.get(config, ())
        records.append(
            {
                "index": index,
                "config": _relative_to_repo(config),
                "config_sha256": _sha256(config),
                "kind": kind,
                "condition_id": condition_id,
                "overrides": list(overrides),
                "command": (
                    _generic_command(config, run_dir, condition_id, overrides)
                    if condition_id is not None
                    else command_for_config(config, run_dir, overrides)
                ),
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
        "aggregation": {"required": job.is_generic_campaign_job, "status": "pending"},
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


def _read_json(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise RuntimeError(f"required artifact is missing: {path}")
    data = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(data, dict):
        raise RuntimeError(f"artifact must be a JSON object: {path}")
    return data


def _aggregate_generic_campaign(job: ParallelJob, manifest: dict[str, Any]) -> Path:
    """Build a canonical campaign view only after every shard completed.

    The condition directories are relative symlinks to child outputs, avoiding a
    second copy of large archives on cs10.  They are deliberately not created
    for partial or failed jobs.
    """

    if manifest["failed_configs"]:
        raise RuntimeError("cannot aggregate a generic campaign with failed shards")
    root = Path(manifest["output_root"])
    config = job.configs[0]
    overrides = list(job.config_overrides.get(config, ()))
    campaign = apply_campaign_cli_overrides(load_yaml(config), overrides)
    all_conditions = {
        item["condition_id"]: item for item in build_campaign_conditions(campaign)
    }
    selected = [all_conditions[condition_id] for condition_id in job.condition_ids]
    base_config_path = Path(str(campaign["base_config"]))
    if not base_config_path.is_absolute():
        base_config_path = REPO_ROOT / base_config_path
    base_cfg = load_yaml(base_config_path)
    base_simulation = SimulationConfig.from_dict(base_cfg).with_overrides(
        campaign.get("base_overrides", {})
    )
    campaign_root = root / "campaign"
    conditions_root = campaign_root / "conditions"
    conditions_root.mkdir(parents=True, exist_ok=False)

    rows: list[dict[str, Any]] = []
    manifests: list[dict[str, Any]] = []
    for record, condition in zip(manifest["configs"], selected, strict=True):
        if record["status"] != "succeeded":
            raise RuntimeError(f"shard did not succeed: {record['index']}")
        child_root = Path(record["output_dir"])
        completion = _read_json(child_root / "campaign_completion.json")
        if completion.get("status") != "completed" or completion.get("exit_code") != 0:
            raise RuntimeError(
                f"shard completion is not successful: {condition['condition_id']}"
            )
        child_manifest = _read_json(child_root / "run_manifest.json")
        child_conditions = child_manifest.get("conditions")
        if not isinstance(child_conditions, list) or len(child_conditions) != 1:
            raise RuntimeError(
                f"shard does not contain exactly one condition: {condition['condition_id']}"
            )
        if child_conditions[0].get("condition_id") != condition["condition_id"]:
            raise RuntimeError(
                f"shard condition ID mismatch: {condition['condition_id']}"
            )
        child_dir = child_root / condition["condition_id"]
        summary = _read_json(child_dir / "run_summary.json")
        if summary.get("execution", {}).get("status") != "completed":
            raise RuntimeError(f"condition is incomplete: {condition['condition_id']}")
        link = conditions_root / condition["condition_id"]
        link.symlink_to(os.path.relpath(child_dir, link.parent))
        cfg = SimulationConfig.from_dict(base_cfg).with_overrides(
            condition["config_overrides"]
        )
        rows.append(_condition_row(cfg, condition, link))
        manifests.append(
            _manifest_condition_record(
                campaign_root,
                condition,
                time_manifest=cfg.time_manifest(),
            )
        )

    summary_path = campaign_root / "summary.csv"
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=_summary_fieldnames(rows))
        writer.writeheader()
        writer.writerows(rows)
    run_manifest = {
        "kind": "generic_multi_run",
        "parallel_aggregate": True,
        "created_at": _now(),
        "campaign_config": _relative_to_repo(config),
        "campaign_overrides": overrides,
        "base_config": str(campaign["base_config"]),
        "source_config_path": str(base_config_path),
        "model_profile": base_simulation.model_profile_manifest(),
        "time": base_simulation.time_manifest(),
        "summary_csv": str(summary_path),
        "output_root": str(campaign_root),
        "axes": campaign_axes_metadata(campaign),
        "condition_order": list(job.condition_ids),
        "conditions": manifests,
        "parallel_job_manifest": str(root / "job_manifest.json"),
    }
    run_manifest.update(base_simulation.implementation_manifest())
    for name, payload in (
        ("run_manifest.json", run_manifest),
        ("manifest.json", run_manifest),
    ):
        (campaign_root / name).write_text(
            json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
    (campaign_root / "campaign_completion.json").write_text(
        json.dumps(
            {
                "status": "completed",
                "expected_condition_count": len(selected),
                "completed_condition_count": len(rows),
                "exit_code": 0,
                "summary_csv": str(summary_path),
                "run_manifest_json": str(campaign_root / "run_manifest.json"),
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    (campaign_root / "run.log").write_text(
        "parallel generic campaign aggregation completed\n", encoding="utf-8"
    )
    return campaign_root


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
    if manifest["status"] == "succeeded" and job.is_generic_campaign_job:
        try:
            campaign_root = _aggregate_generic_campaign(job, manifest)
        except Exception as exc:
            manifest["status"] = "failed"
            manifest["aggregation"] = {
                "required": True,
                "status": "failed",
                "error": str(exc),
            }
        else:
            manifest["aggregation"] = {
                "required": True,
                "status": "completed",
                "campaign_root": str(campaign_root),
            }
    elif job.is_generic_campaign_job:
        manifest["aggregation"] = {
            "required": True,
            "status": "not_created_due_to_shard_failure",
        }
    with (root / "job.log").open("a", encoding="utf-8") as handle:
        handle.write(f"parallel job {manifest['status']}\n")
    _write_manifest(root, manifest)
    return manifest
