from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import pytest

from sim_swim.analysis.parallel_job import (
    ParallelJob,
    build_plan,
    load_parallel_job,
    resolve_execution,
    run_parallel_job,
)


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "conf/phase2_parallel/example_stage_a_validation/job.yaml"
SWEEP_A = ROOT / "conf/phase2_sweeps/2015_stage_a_motor_off.yaml"
SWEEP_B = ROOT / "conf/phase2_sweeps/2015_stage_a_motor_on.yaml"
SHAPE_SWEEP = ROOT / "conf/phase2_sweeps/shape_stability_grid.yaml"


def _job() -> ParallelJob:
    return ParallelJob(
        schema_version=1,
        job_id="test_job",
        job_name="test_job",
        config_path=EXAMPLE.resolve(),
        configs=(SWEEP_A.resolve(), SWEEP_B.resolve()),
        max_workers="auto",
        worker_policy="cs10_qualified",
    )


def _load_script(name: str):
    path = ROOT / "scripts/01_simulate_swimming/run_parallel.py"
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def test_load_example_job_reuses_existing_sweep_profiles() -> None:
    job = load_parallel_job(EXAMPLE)

    assert job.job_name == "example_stage_a_validation"
    assert job.configs == (SWEEP_A.resolve(), SWEEP_B.resolve())
    assert job.max_workers == "auto"


def test_worker_validation_rejects_invalid_override() -> None:
    with pytest.raises(ValueError, match="positive integer"):
        resolve_execution(_job(), 0)


def test_auto_worker_policies_and_thread_environment(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    job = _job()
    resolved = resolve_execution(job, None)
    assert resolved.max_workers == 2
    assert resolved.thread_environment == {
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }

    host_job = ParallelJob(**{**job.__dict__, "worker_policy": "host_cpu"})
    monkeypatch.setattr("sim_swim.analysis.parallel_job.os.cpu_count", lambda: 1)
    assert resolve_execution(host_job, "auto").max_workers == 1
    assert resolve_execution(host_job, 2).thread_environment == {}


def test_build_plan_uses_stage_a_output_base_dir(tmp_path: Path) -> None:
    plan = build_plan(_job(), resolve_execution(_job(), None), tmp_path / "root")

    command = plan["configs"][0]["command"]
    assert "config=conf/phase2_sweeps/2015_stage_a_motor_off.yaml" in command
    assert any(item.startswith("output_base_dir=") for item in command)
    assert plan["configs"][0]["output_dir"].endswith("001_2015_stage_a_motor_off/run")


def test_build_plan_uses_output_dir_for_non_stage_sweep(tmp_path: Path) -> None:
    job = ParallelJob(
        schema_version=1,
        job_id="shape",
        job_name="shape",
        config_path=EXAMPLE.resolve(),
        configs=(SHAPE_SWEEP.resolve(),),
        max_workers=1,
        worker_policy="host_cpu",
    )
    command = build_plan(job, resolve_execution(job, None), tmp_path)["configs"][0][
        "command"
    ]

    assert any(item.startswith("output_dir=") for item in command)
    assert not any(item.startswith("output_base_dir=") for item in command)


class _FakeProcess:
    def __init__(self, exit_code: int) -> None:
        self.exit_code = exit_code

    def poll(self) -> int:
        return self.exit_code


def test_partial_failure_keeps_other_config_records(tmp_path: Path) -> None:
    outcomes = iter([0, 7])

    def fake_popen(*args: object, **kwargs: object) -> _FakeProcess:
        assert kwargs["cwd"] == ROOT
        assert kwargs["env"]["OMP_NUM_THREADS"] == "1"
        return _FakeProcess(next(outcomes))

    manifest = run_parallel_job(
        _job(),
        resolve_execution(_job(), None),
        output_base_dir=tmp_path,
        popen=fake_popen,
        poll_interval_s=0,
    )

    assert manifest["status"] == "failed"
    assert manifest["failed_configs"] == [2]
    assert [item["status"] for item in manifest["configs"]] == ["succeeded", "failed"]
    assert all("wall_time_s" in item for item in manifest["configs"])
    saved = json.loads(
        (Path(manifest["output_root"]) / "job_manifest.json").read_text()
    )
    assert saved["completion_order"] == [1, 2]
    assert (Path(manifest["configs"][0]["stdout_log"])).is_file()
    assert (Path(manifest["configs"][1]["stderr_log"])).is_file()


def test_failed_process_start_is_recorded_and_next_config_runs(tmp_path: Path) -> None:
    outcomes = iter([OSError("unavailable"), _FakeProcess(0)])

    def fake_popen(*args: object, **kwargs: object) -> _FakeProcess:
        outcome = next(outcomes)
        if isinstance(outcome, OSError):
            raise outcome
        return outcome

    manifest = run_parallel_job(
        _job(),
        resolve_execution(_job(), 1),
        output_base_dir=tmp_path,
        popen=fake_popen,
        poll_interval_s=0,
    )

    assert manifest["failed_configs"] == [1]
    assert [item["status"] for item in manifest["configs"]] == [
        "failed_to_start",
        "succeeded",
    ]


def test_dry_run_prints_plan_without_creating_output(
    capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    launcher = _load_script("phase2_run_parallel_dry_run")
    expected_root = tmp_path / "planned"
    monkeypatch.setattr(launcher, "job_output_root", lambda job: expected_root)

    assert launcher.main([f"config={EXAMPLE.relative_to(ROOT)}", "dry_run=true"]) == 0

    plan = json.loads(capsys.readouterr().out)
    assert plan["status"] == "planned"
    assert plan["output_root"] == str(expected_root)
    assert plan["execution"]["thread_environment_overrides"] == {
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    }
    assert not expected_root.exists()
