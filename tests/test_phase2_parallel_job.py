from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys

import pytest

from sim_swim.analysis import parallel_job
from sim_swim.analysis.parallel_job import (
    ParallelJob,
    build_plan,
    load_parallel_job,
    resolve_execution,
    run_parallel_job,
)


ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "conf/phase2_parallel/example_stage_a_validation/job.yaml"
ISSUE210_2010 = ROOT / "conf/phase2_parallel/issue210_2010_project/job.yaml"
ISSUE203 = ROOT / "conf/phase2_parallel/issue203_uniform_torque_profile/job.yaml"
ISSUE203_QUALIFICATION = (
    ROOT / "conf/phase2_parallel/issue203_uniform_torque_profile/qualification_job.yaml"
)
ISSUE203_DT_CONTACT = (
    ROOT / "conf/phase2_parallel/issue203_torque_profile_dt_contact/job.yaml"
)
ISSUE215 = ROOT / "conf/phase2_parallel/issue215_5s_axis_convergence/job.yaml"
ISSUE215_QUALIFICATION = (
    ROOT / "conf/phase2_parallel/issue215_5s_axis_convergence/qualification_job.yaml"
)
SWEEP_A = ROOT / "conf/phase2_sweeps/2015_stage_a_motor_off.yaml"
SWEEP_B = ROOT / "conf/phase2_sweeps/2015_stage_a_motor_on.yaml"
SHAPE_SWEEP = ROOT / "conf/phase2_sweeps/shape_stability_grid.yaml"
TORQUE_DISTRIBUTION_SWEEP = ROOT / "conf/phase2_sweeps/torque_distribution_grid.yaml"


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


def test_issue210_job_uses_existing_2010_project_profiles() -> None:
    job = load_parallel_job(ISSUE210_2010)

    assert [config.name for config in job.configs] == [
        "shape_stability_grid.yaml",
        "torque_distribution_grid.yaml",
    ]
    assert job.worker_policy == "cs10_qualified"
    assert job.config_overrides == {
        SHAPE_SWEEP.resolve(): (
            "sample_limit=1",
            "attach_seed=0",
            "phase_seed=0",
            "duration_s=0.001",
        ),
        TORQUE_DISTRIBUTION_SWEEP.resolve(): (
            "sample_limit=1",
            "torque_nm=2.5e-20",
            "duration_s=0.001",
        ),
    }


def test_issue203_generic_job_expands_27_independent_conditions() -> None:
    job = load_parallel_job(ISSUE203)
    execution = resolve_execution(job, None)
    plan = build_plan(job, execution, ROOT / ".tmp_issue203_plan")

    assert job.is_generic_campaign_job
    assert job.task_count == 27
    assert execution.max_workers == 8
    assert len(plan["configs"]) == 27
    assert plan["configs"][0]["condition_id"] == "as000__ps000__nf01"
    assert plan["configs"][-1]["condition_id"] == "as002__ps002__nf03"
    command = plan["configs"][0]["command"]
    assert "scripts/01_simulate_swimming/run_multi_run.py" in command
    assert "output.timestamp_subdir=false" in command
    assert "sweep.include_condition_ids=[as000__ps000__nf01]" in command


def test_issue203_qualification_job_preserves_27_shards_and_duration_override() -> None:
    job = load_parallel_job(ISSUE203_QUALIFICATION)
    plan = build_plan(job, resolve_execution(job, None), ROOT / ".tmp_issue203_plan")

    assert job.task_count == 27
    assert all(
        record["overrides"] == ["time.duration_s=0.001"] for record in plan["configs"]
    )
    assert all(
        "time.duration_s=0.001" in record["command"] for record in plan["configs"]
    )


def test_issue203_dt_contact_job_expands_only_the_12_new_shards() -> None:
    job = load_parallel_job(ISSUE203_DT_CONTACT)
    plan = build_plan(job, resolve_execution(job, None), ROOT / ".tmp_issue203_plan")

    assert job.task_count == 12
    assert len(plan["configs"]) == 12
    assert all("nf03" in str(record["condition_id"]) for record in plan["configs"])
    assert all(
        "dt1e-3" not in str(record["condition_id"]) for record in plan["configs"]
    )


def test_issue215_job_expands_36_issue204_matching_conditions() -> None:
    job = load_parallel_job(ISSUE215)
    execution = resolve_execution(job, None)
    plan = build_plan(job, execution, ROOT / ".tmp_issue215_plan")

    assert job.is_generic_campaign_job
    assert job.task_count == 36
    assert execution.max_workers == 8
    assert len(plan["configs"]) == 36
    assert plan["configs"][0]["condition_id"] == "as000__ps000__nf01"
    assert plan["configs"][-1]["condition_id"] == "as002__ps002__nf04"


def test_issue215_qualification_preserves_36_shards_and_duration_override() -> None:
    job = load_parallel_job(ISSUE215_QUALIFICATION)
    plan = build_plan(job, resolve_execution(job, None), ROOT / ".tmp_issue215_plan")

    assert job.task_count == 36
    assert all(
        record["overrides"] == ["time.duration_s=0.001"] for record in plan["configs"]
    )


def test_generic_aggregate_requires_all_shards_and_creates_canonical_view(
    tmp_path: Path,
) -> None:
    job = load_parallel_job(ISSUE203)
    execution = resolve_execution(job, None)
    root = tmp_path / "parallel_job"
    manifest = build_plan(job, execution, root)
    root.mkdir()
    manifest["output_root"] = str(root)
    for record in manifest["configs"]:
        child_root = Path(record["output_dir"])
        condition_dir = child_root / str(record["condition_id"])
        condition_dir.mkdir(parents=True)
        (child_root / "campaign_completion.json").write_text(
            json.dumps({"status": "completed", "exit_code": 0}), encoding="utf-8"
        )
        (child_root / "run_manifest.json").write_text(
            json.dumps({"conditions": [{"condition_id": record["condition_id"]}]}),
            encoding="utf-8",
        )
        (condition_dir / "run_summary.json").write_text(
            json.dumps(
                {
                    "execution": {"status": "completed"},
                    "gates": {},
                    "all_step_metrics": {},
                }
            ),
            encoding="utf-8",
        )
        record["status"] = "succeeded"

    campaign = parallel_job._aggregate_generic_campaign(job, manifest)

    assert (campaign / "campaign_completion.json").is_file()
    assert (campaign / "summary.csv").is_file()
    first_link = campaign / "conditions/as000__ps000__nf01"
    assert first_link.is_symlink()
    aggregate_manifest = json.loads((campaign / "run_manifest.json").read_text())
    assert aggregate_manifest["conditions"][0]["output_dir"] == str(first_link)
    assert Path(aggregate_manifest["conditions"][0]["output_dir"]).is_dir()

    manifest["failed_configs"] = [1]
    with pytest.raises(RuntimeError, match="failed shards"):
        parallel_job._aggregate_generic_campaign(job, manifest)


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        ("      - duration_s\n", "config overrides must use KEY=VALUE"),
        ("      - output_dir=bad\n", "must not set launcher-managed output paths"),
        ("      - output-dir=bad\n", "must not set launcher-managed output paths"),
        ("      - output-base-dir=bad\n", "must not set launcher-managed output paths"),
        ("      - duration_s=0.001\n      - duration_s=0.002\n", "must not repeat key"),
        ("      - duration_s=0.001\n      - duration-s=0.002\n", "must not repeat key"),
    ],
)
def test_job_rejects_invalid_config_overrides(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, overrides: str, message: str
) -> None:
    sweep_dir = tmp_path / "conf/phase2_sweeps"
    sweep_dir.mkdir(parents=True)
    (sweep_dir / "sweep.yaml").write_text(
        "kind: shape_stability_grid\nmetadata:\n  role: sweep\nargs: {}\n",
        encoding="utf-8",
    )
    job_path = tmp_path / "conf/phase2_parallel/test/job.yaml"
    job_path.parent.mkdir(parents=True)
    job_path.write_text(
        "schema_version: 1\njob_id: test\nconfigs:\n"
        "  - path: conf/phase2_sweeps/sweep.yaml\n    overrides:\n"
        f"{overrides}execution:\n  max_workers: 1\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(parallel_job, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(parallel_job, "SWEEP_DIRECTORY", sweep_dir)

    with pytest.raises(ValueError, match=message):
        load_parallel_job(job_path)


@pytest.mark.parametrize(
    ("job_extra", "execution_extra", "message"),
    [
        ("unexpected: true\n", "", "job config contains unsupported keys: unexpected"),
        (
            "",
            "  worker_polciy: cs10_qualified\n",
            "execution contains unsupported keys: worker_polciy",
        ),
    ],
)
def test_job_rejects_unknown_schema_keys(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    job_extra: str,
    execution_extra: str,
    message: str,
) -> None:
    sweep_dir = tmp_path / "conf/phase2_sweeps"
    sweep_dir.mkdir(parents=True)
    (sweep_dir / "sweep.yaml").write_text(
        "kind: shape_stability_grid\nmetadata:\n  role: sweep\nargs: {}\n",
        encoding="utf-8",
    )
    job_path = tmp_path / "conf/phase2_parallel/test/job.yaml"
    job_path.parent.mkdir(parents=True)
    job_path.write_text(
        "schema_version: 1\n"
        "job_id: test\n"
        "configs:\n"
        "  - conf/phase2_sweeps/sweep.yaml\n"
        f"{job_extra}"
        "execution:\n"
        "  max_workers: auto\n"
        f"{execution_extra}",
        encoding="utf-8",
    )
    monkeypatch.setattr(parallel_job, "REPO_ROOT", tmp_path)
    monkeypatch.setattr(parallel_job, "SWEEP_DIRECTORY", sweep_dir)

    with pytest.raises(ValueError, match=message):
        load_parallel_job(job_path)


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

    oversized_cs10_job = ParallelJob(
        **{**job.__dict__, "configs": (SWEEP_A.resolve(),) * 9}
    )
    assert resolve_execution(oversized_cs10_job, 20).max_workers == 8


def test_job_output_root_uses_jst_date_and_time_namespace(tmp_path: Path) -> None:
    root = parallel_job.job_output_root(_job(), output_base_dir=tmp_path)

    assert root.parent.parent.parent.parent == tmp_path
    assert root.parent.parent.name.isdigit() and len(root.parent.parent.name) == 6
    assert root.parent.parent.parent.name.count("-") == 2
    assert root.parent.name == "parallel"
    assert root.name.startswith("test_job__")


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


def test_build_plan_records_config_overrides(tmp_path: Path) -> None:
    job = ParallelJob(
        schema_version=1,
        job_id="shape",
        job_name="shape",
        config_path=EXAMPLE.resolve(),
        configs=(SHAPE_SWEEP.resolve(),),
        max_workers=1,
        worker_policy="host_cpu",
        config_overrides={
            SHAPE_SWEEP.resolve(): ("sample_limit=1", "duration_s=0.001")
        },
    )

    record = build_plan(job, resolve_execution(job, None), tmp_path)["configs"][0]

    assert record["overrides"] == ["sample_limit=1", "duration_s=0.001"]
    assert record["command"][-2:] == ["sample_limit=1", "duration_s=0.001"]


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


def test_requested_output_root_is_used_without_timestamp_generation(
    tmp_path: Path,
) -> None:
    root = tmp_path / "fixed_parallel_root"

    manifest = run_parallel_job(
        _job(),
        resolve_execution(_job(), 2),
        output_root=root,
        popen=lambda *args, **kwargs: _FakeProcess(0),
        poll_interval_s=0,
    )

    assert Path(manifest["output_root"]) == root
    assert (root / "job_manifest.json").is_file()


def test_requested_output_root_rejects_output_base_dir(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="output_base_dir or output_root"):
        run_parallel_job(
            _job(),
            resolve_execution(_job(), 1),
            output_base_dir=tmp_path,
            output_root=tmp_path / "fixed",
        )


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
