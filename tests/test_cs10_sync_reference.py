from __future__ import annotations

import hashlib
import importlib.util
from pathlib import Path
import shutil
import sys
from types import SimpleNamespace

import pytest


ROOT = Path(__file__).resolve().parents[1]
SYNC_SCRIPT = ROOT / "scripts/cs10/sync_reference_from_cs10.py"


def _module():
    spec = importlib.util.spec_from_file_location(
        "cs10_sync_reference_test", SYNC_SCRIPT
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _campaign(tmp_path: Path) -> tuple[Path, Path]:
    campaign = tmp_path / "remote/campaign"
    child = tmp_path / "remote/children/nf01__slots0"
    child.mkdir(parents=True)
    (child / "state_archive.npz").write_bytes(b"completed archive")
    (child / "run_summary.json").write_bytes(b"summary")
    (child / "performance.json").write_bytes(b"performance")
    (child / "run.log").write_bytes(b"operational log")
    conditions = campaign / "conditions"
    conditions.mkdir(parents=True)
    (conditions / "nf01__slots0").symlink_to(child, target_is_directory=True)
    for name in (
        "run_manifest.json",
        "manifest.json",
        "summary.csv",
        "campaign_completion.json",
    ):
        (campaign / name).write_bytes(name.encode())
    return campaign, child


def test_parallel_campaign_dry_run_selects_only_live_artifacts(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    module = _module()
    commands: list[list[str]] = []
    monkeypatch.setattr(module, "_run", lambda args, *, dry_run: commands.append(args))
    destination = tmp_path / "local"
    module.sync(
        "cs10",
        "/remote/campaign",
        destination,
        layout="parallel-campaign",
        dry_run=True,
    )
    assert not destination.exists()
    assert commands[0][:4] == [
        "rsync",
        "-aL",
        "--exclude=run.log",
        "--exclude=render.log",
    ]
    assert commands[0][4] == "cs10:/remote/campaign/conditions/"
    assert [Path(command[1]).name for command in commands[1:]] == list(
        module.PARALLEL_CAMPAIGN_ROOT_FILES
    )
    assert not any(
        "/analysis" in item or "reference_manifest" in item
        for command in commands
        for item in command
    )


def test_reference_layout_retains_existing_transfer_list(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    module = _module()
    commands: list[list[str]] = []
    monkeypatch.setattr(module, "_run", lambda args, *, dry_run: commands.append(args))
    module.sync("cs10", "/remote/reference", tmp_path / "local", dry_run=True)
    assert commands[0][:3] == ["scp", "-r", "cs10:/remote/reference/conditions"]
    assert commands[1][:3] == ["scp", "-r", "cs10:/remote/reference/analysis"]
    assert any(
        "reference_manifest.json" in item for command in commands for item in command
    )


def test_remote_hashes_follow_parallel_symlinks_and_exclude_logs(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    module = _module()
    commands: list[list[str]] = []

    def fake_run(args: list[str], **_kwargs: object) -> SimpleNamespace:
        commands.append(args)
        return SimpleNamespace(
            stdout="abc123  conditions/nf01__slots0/state_archive.npz\n"
        )

    monkeypatch.setattr(module.subprocess, "run", fake_run)
    result = module._remote_hashes(
        "cs10", "/remote/campaign", layout="parallel-campaign"
    )
    assert result == {"conditions/nf01__slots0/state_archive.npz": "abc123"}
    command = commands[0][2]
    assert (
        "find -L conditions run_manifest.json manifest.json summary.csv campaign_completion.json"
        in command
    )
    assert "! -name run.log ! -name render.log" in command


def test_parallel_campaign_sync_dereferences_and_verifies_hashes(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    module = _module()
    campaign, child = _campaign(tmp_path)
    local = tmp_path / "local"

    def fake_transfer(args: list[str], *, dry_run: bool) -> None:
        assert not dry_run
        if args[0] == "rsync":
            shutil.copytree(
                campaign / "conditions",
                local / "conditions",
                symlinks=False,
                dirs_exist_ok=True,
                ignore=shutil.ignore_patterns("run.log", "render.log"),
            )
        else:
            source = campaign / Path(args[1]).name
            if not source.is_file():
                raise FileNotFoundError(source)
            shutil.copy2(source, local)

    def remote_hashes(_host: str, _remote_dir: str, *, layout: str) -> dict[str, str]:
        assert layout == "parallel-campaign"
        files = {
            f"conditions/nf01__slots0/{name}": child / name
            for name in ("state_archive.npz", "run_summary.json", "performance.json")
        }
        files.update(
            {name: campaign / name for name in module.PARALLEL_CAMPAIGN_ROOT_FILES}
        )
        return {
            name: hashlib.sha256(path.read_bytes()).hexdigest()
            for name, path in files.items()
        }

    monkeypatch.setattr(module, "_run", fake_transfer)
    monkeypatch.setattr(module, "_remote_hashes", remote_hashes)
    assert (
        module.sync("cs10", str(campaign), local, layout="parallel-campaign") == local
    )
    assert (
        local / "conditions/nf01__slots0/state_archive.npz"
    ).read_bytes() == b"completed archive"
    assert not (local / "conditions/nf01__slots0").is_symlink()
    assert not (local / "conditions/nf01__slots0/run.log").exists()


def test_parallel_campaign_sync_rejects_missing_artifact_and_hash_mismatch(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    module = _module()
    campaign, _child = _campaign(tmp_path)

    def fake_transfer(args: list[str], *, dry_run: bool) -> None:
        assert not dry_run
        destination = Path(args[-1])
        if args[0] == "rsync":
            shutil.copytree(
                campaign / "conditions", destination, symlinks=False, dirs_exist_ok=True
            )
        else:
            source = campaign / Path(args[1]).name
            if not source.is_file():
                raise FileNotFoundError(source)
            shutil.copy2(source, destination)

    monkeypatch.setattr(module, "_run", fake_transfer)
    (campaign / "summary.csv").unlink()
    with pytest.raises(FileNotFoundError, match="summary.csv"):
        module.sync(
            "cs10", str(campaign), tmp_path / "missing", layout="parallel-campaign"
        )

    (campaign / "summary.csv").write_bytes(b"summary.csv")
    monkeypatch.setattr(
        module,
        "_remote_hashes",
        lambda *_args, **_kwargs: {"manifest.json": "0" * 64},
    )
    with pytest.raises(RuntimeError, match="checksum mismatch"):
        module.sync(
            "cs10", str(campaign), tmp_path / "mismatch", layout="parallel-campaign"
        )
