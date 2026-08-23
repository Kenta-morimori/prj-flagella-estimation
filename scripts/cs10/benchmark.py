#!/usr/bin/env python3
"""Benchmark fixed independent Phase 2 workloads on cs10.

This is deliberately not a general parallel-job launcher. Issue #209 owns that
interface; this runner only establishes the safe worker count for it.
"""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from datetime import datetime
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from typing import Any, Iterable
from zoneinfo import ZoneInfo

import yaml


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = REPOSITORY_ROOT / "src"
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))


DEFAULT_WORKERS = (1, 2, 4, 6, 8, 10)
THREAD_ENV = {
    "OMP_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
}
TIME_RSS_PATTERN = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")
TIME_CPU_PATTERN = re.compile(r"Percent of CPU this job got:\s*(\d+)%")


@dataclass(frozen=True)
class JobResult:
    worker_count: int
    repetition: int
    job_index: int
    returncode: int
    elapsed_s: float
    peak_rss_kib: int | None
    cpu_percent: int | None
    output_dir: str
    stdout_log: str
    time_log: str


def _now_output_root() -> Path:
    now = datetime.now(ZoneInfo("Asia/Tokyo"))
    return (
        Path("outputs")
        / now.strftime("%Y-%m-%d")
        / now.strftime("%H%M%S")
        / "cs10_qualification"
    )


def _parse_worker_counts(value: str) -> tuple[int, ...]:
    values = tuple(int(item) for item in value.split(",") if item)
    if not values or any(item < 1 for item in values):
        raise argparse.ArgumentTypeError("worker counts must be positive integers")
    return values


def _median(values: Iterable[float]) -> float | None:
    ordered = sorted(values)
    if not ordered:
        return None
    mid = len(ordered) // 2
    return ordered[mid] if len(ordered) % 2 else (ordered[mid - 1] + ordered[mid]) / 2


def recommend_worker_count(
    rows: list[dict[str, Any]], *, reserve_cores: int = 2
) -> int | None:
    """Return the fastest non-failing setting that leaves physical cores reserved."""
    limit = 10 - reserve_cores
    candidates = [
        row
        for row in rows
        if row["worker_count"] <= limit
        and row["successes"] == row["worker_count"] * row["repetitions"]
    ]
    if not candidates:
        return None
    return max(
        candidates,
        key=lambda row: (row["median_throughput_jobs_per_s"], row["worker_count"]),
    )["worker_count"]


def _read_time_metrics(path: Path) -> tuple[int | None, int | None]:
    if not path.is_file():
        return None, None
    text = path.read_text(encoding="utf-8", errors="replace")
    rss = TIME_RSS_PATTERN.search(text)
    cpu = TIME_CPU_PATTERN.search(text)
    return (int(rss.group(1)) if rss else None, int(cpu.group(1)) if cpu else None)


def workload_command(*, output_dir: Path, duration_s: float) -> list[str]:
    """Build the fixed qualification workload command without launching it."""
    return [
        sys.executable,
        "scripts/01_simulate_swimming/run_sweep.py",
        "config=conf/phase2_sweeps/shape_stability_grid.yaml",
        "sample_limit=1",
        "attach_seed=0",
        "phase_seed=0",
        f"duration_s={duration_s}",
        f"output_dir={output_dir / 'campaign'}",
        "overwrite=true",
    ]


def workload_metadata(repository_root: Path, duration_s: float) -> dict[str, Any]:
    """Resolve the effective integration timing of the fixed workload."""
    from sim_swim.sim.params import SimulationConfig

    profile_path = repository_root / "conf/phase2_sweeps/shape_stability_grid.yaml"
    profile = yaml.safe_load(profile_path.read_text(encoding="utf-8")) or {}
    args = profile["args"]
    base_config_path = repository_root / args["config"]
    base_config = yaml.safe_load(base_config_path.read_text(encoding="utf-8")) or {}
    torque_nm = float(args["torque_nm"])
    config = SimulationConfig.from_dict(base_config).with_overrides(
        {
            "flagella": {"n_flagella": int(args["n_flagella"])},
            "motor": {
                "torque_Nm": torque_nm,
                "reference_torque_Nm": torque_nm,
            },
            "time": {
                "duration_s": duration_s,
                "dt_star": float(args["dt_star"]),
            },
        }
    )
    return {
        "config": str(profile_path.relative_to(repository_root)),
        "base_config": str(base_config_path.relative_to(repository_root)),
        "sample_limit": 1,
        "n_flagella": config.flagella.n_flagella,
        "duration_s": config.time.duration_s,
        "dt_star": config.dt_star,
        "tau_s": config.tau_s,
        "dt_s": config.dt_s,
        "total_steps": config.total_steps,
        "state_archive": True,
    }


def _run_job(
    *,
    repository_root: Path,
    output_dir: Path,
    worker_count: int,
    repetition: int,
    job_index: int,
    thread_limited: bool,
    duration_s: float,
) -> JobResult:
    output_dir.mkdir(parents=True, exist_ok=False)
    stdout_path = output_dir / "stdout.log"
    time_path = output_dir / "time.txt"
    command = workload_command(output_dir=output_dir, duration_s=duration_s)
    time_binary = shutil.which("/usr/bin/time") or shutil.which("time")
    if time_binary:
        command = [time_binary, "-v", "-o", str(time_path), *command]
    environment = os.environ.copy()
    if thread_limited:
        environment.update(THREAD_ENV)
    else:
        for name in THREAD_ENV:
            environment.pop(name, None)
    start = time.monotonic()
    with stdout_path.open("w", encoding="utf-8") as handle:
        completed = subprocess.run(
            command,
            cwd=repository_root,
            env=environment,
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    elapsed_s = time.monotonic() - start
    peak_rss_kib, cpu_percent = _read_time_metrics(time_path)
    return JobResult(
        worker_count=worker_count,
        repetition=repetition,
        job_index=job_index,
        returncode=completed.returncode,
        elapsed_s=elapsed_s,
        peak_rss_kib=peak_rss_kib,
        cpu_percent=cpu_percent,
        output_dir=str(output_dir),
        stdout_log=str(stdout_path),
        time_log=str(time_path),
    )


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fieldnames = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _recommendations(summary_rows: list[dict[str, Any]]) -> dict[str, int | None]:
    return {
        "threads_unset": recommend_worker_count(
            [row for row in summary_rows if row["thread_limited"] == "0"]
        ),
        "threads_1": recommend_worker_count(
            [row for row in summary_rows if row["thread_limited"] == "1"]
        ),
    }


def summarize_existing(output_dir: Path) -> dict[str, int | None]:
    """Recompute recommendations from an existing screen without rerunning jobs."""
    with (output_dir / "summary.csv").open(encoding="utf-8", newline="") as handle:
        rows = [
            {
                **row,
                "worker_count": int(row["worker_count"]),
                "repetitions": int(row["repetitions"]),
                "successes": int(row["successes"]),
                "median_throughput_jobs_per_s": float(
                    row["median_throughput_jobs_per_s"]
                ),
            }
            for row in csv.DictReader(handle)
        ]
    recommendations = _recommendations(rows)
    manifest_path = output_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["recommendations"] = recommendations
    manifest["recommendations_recomputed_at"] = datetime.now(
        ZoneInfo("Asia/Tokyo")
    ).isoformat()
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return recommendations


def run_benchmark(args: argparse.Namespace) -> Path:
    root = args.output_dir
    root.mkdir(parents=True, exist_ok=False)
    records: list[JobResult] = []
    for thread_limited in (False, True):
        for worker_count in args.worker_counts:
            for repetition in range(1, args.repetitions + 1):
                started = time.monotonic()
                with ThreadPoolExecutor(max_workers=worker_count) as executor:
                    futures = [
                        executor.submit(
                            _run_job,
                            repository_root=args.repository_root,
                            output_dir=root
                            / ("threads_1" if thread_limited else "threads_unset")
                            / f"workers_{worker_count}"
                            / f"repeat_{repetition}"
                            / f"job_{index:02d}",
                            worker_count=worker_count,
                            repetition=repetition,
                            job_index=index,
                            thread_limited=thread_limited,
                            duration_s=args.duration_s,
                        )
                        for index in range(1, worker_count + 1)
                    ]
                    records.extend(future.result() for future in futures)
                print(
                    f"threads_limited={thread_limited} workers={worker_count} repetition={repetition} elapsed_s={time.monotonic() - started:.3f}",
                    flush=True,
                )

    record_rows = [
        record.__dict__
        | {"thread_limited": "1" if "threads_1" in record.output_dir else "0"}
        for record in records
    ]
    _write_csv(root / "job_results.csv", record_rows)
    summary_rows: list[dict[str, Any]] = []
    for thread_limited in (False, True):
        for worker_count in args.worker_counts:
            subset = [
                record
                for record in records
                if record.worker_count == worker_count
                and ("threads_1" in record.output_dir) == thread_limited
            ]
            repetition_throughputs = []
            for repetition in range(1, args.repetitions + 1):
                run_rows = [
                    record for record in subset if record.repetition == repetition
                ]
                elapsed = max((record.elapsed_s for record in run_rows), default=0.0)
                if elapsed:
                    repetition_throughputs.append(worker_count / elapsed)
            summary_rows.append(
                {
                    "thread_limited": "1" if thread_limited else "0",
                    "worker_count": worker_count,
                    "repetitions": args.repetitions,
                    "successes": sum(record.returncode == 0 for record in subset),
                    "median_throughput_jobs_per_s": _median(repetition_throughputs)
                    or 0.0,
                    "median_peak_rss_kib": _median(
                        [
                            float(record.peak_rss_kib)
                            for record in subset
                            if record.peak_rss_kib is not None
                        ]
                    )
                    or 0.0,
                    "median_cpu_percent": _median(
                        [
                            float(record.cpu_percent)
                            for record in subset
                            if record.cpu_percent is not None
                        ]
                    )
                    or 0.0,
                }
            )
    _write_csv(root / "summary.csv", summary_rows)
    recommendations = _recommendations(summary_rows)
    manifest = {
        "pipeline_name": "cs10_worker_benchmark",
        "created_at": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(),
        "workload": workload_metadata(args.repository_root, args.duration_s),
        "worker_counts": list(args.worker_counts),
        "repetitions": args.repetitions,
        "thread_limited_environment": THREAD_ENV,
        "recommendations": recommendations,
        "outputs": {
            "job_results_csv": str(root / "job_results.csv"),
            "summary_csv": str(root / "summary.csv"),
        },
    }
    (root / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (root / "run.log").write_text(
        f"pipeline_name=cs10_worker_benchmark\nsummary={root / 'summary.csv'}\n",
        encoding="utf-8",
    )
    return root


def main(argv: list[str] | None = None) -> Path:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repository-root", type=Path, default=Path.cwd())
    parser.add_argument("--output-dir", type=Path, default=_now_output_root())
    parser.add_argument(
        "--worker-counts", type=_parse_worker_counts, default=DEFAULT_WORKERS
    )
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--duration-s", type=float, default=0.05)
    parser.add_argument(
        "--summarize-existing",
        type=Path,
        metavar="OUTPUT_DIR",
        help="Recompute recommendations from an existing benchmark output.",
    )
    args = parser.parse_args(argv)
    if args.summarize_existing is not None:
        recommendations = summarize_existing(args.summarize_existing)
        print(json.dumps(recommendations, ensure_ascii=False))
        return args.summarize_existing
    if args.repetitions < 1:
        parser.error("--repetitions must be positive")
    if args.duration_s <= 0:
        parser.error("--duration-s must be positive")
    root = run_benchmark(args)
    print(root)
    return root


if __name__ == "__main__":
    main()
