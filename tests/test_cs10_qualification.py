from __future__ import annotations

import importlib.util
from pathlib import Path
import sys


def _load_script(name: str, relative_path: str):
    path = Path(__file__).resolve().parents[1] / relative_path
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def test_runtime_probe_writes_manifest_and_log(tmp_path: Path) -> None:
    probe = _load_script("cs10_probe", "scripts/cs10/probe_runtime.py")

    output_dir = tmp_path / "probe"
    result = probe.main(["--output-dir", str(output_dir)])

    assert result == output_dir / "manifest.json"
    assert result.is_file()
    assert (output_dir / "run.log").is_file()
    manifest = __import__("json").loads(result.read_text(encoding="utf-8"))
    assert manifest["pipeline_name"] == "cs10_runtime_probe"
    assert manifest["imports"]["sim_swim"]
    assert manifest["packages"]["numpy"]


def test_worker_recommendation_reserves_two_physical_cores() -> None:
    benchmark = _load_script("cs10_benchmark", "scripts/cs10/benchmark.py")
    rows = [
        {
            "worker_count": 4,
            "successes": 5,
            "repetitions": 5,
            "median_throughput_jobs_per_s": 2.0,
        },
        {
            "worker_count": 8,
            "successes": 5,
            "repetitions": 5,
            "median_throughput_jobs_per_s": 3.0,
        },
        {
            "worker_count": 10,
            "successes": 5,
            "repetitions": 5,
            "median_throughput_jobs_per_s": 4.0,
        },
    ]

    assert benchmark.recommend_worker_count(rows) == 8


def test_worker_recommendation_rejects_failures() -> None:
    benchmark = _load_script("cs10_benchmark_fail", "scripts/cs10/benchmark.py")
    rows = [
        {
            "worker_count": 4,
            "successes": 5,
            "repetitions": 5,
            "median_throughput_jobs_per_s": 2.0,
        },
        {
            "worker_count": 8,
            "successes": 4,
            "repetitions": 5,
            "median_throughput_jobs_per_s": 4.0,
        },
    ]

    assert benchmark.recommend_worker_count(rows) == 4


def test_workload_command_records_requested_duration(tmp_path: Path) -> None:
    benchmark = _load_script("cs10_benchmark_duration", "scripts/cs10/benchmark.py")

    command = benchmark.workload_command(output_dir=tmp_path / "job", duration_s=0.001)

    assert "duration_s=0.001" in command
    assert f"output_dir={tmp_path / 'job' / 'campaign'}" in command
