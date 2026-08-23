#!/usr/bin/env python3
"""Record reproducible cs10 system and Python runtime qualification evidence."""

from __future__ import annotations

import argparse
from datetime import datetime
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
from typing import Any
from zoneinfo import ZoneInfo


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = REPOSITORY_ROOT / "src"
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))


PACKAGES = (
    "matplotlib",
    "numpy",
    "opencv-python",
    "pandas",
    "pydantic",
    "PyYAML",
    "rich",
    "scikit-image",
    "scipy",
    "tqdm",
    "typer",
    "pytest",
)
THREAD_VARIABLES = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)


def _command(command: list[str]) -> dict[str, Any]:
    executable = shutil.which(command[0])
    if executable is None:
        return {"available": False, "command": command}
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    return {
        "available": True,
        "command": command,
        "returncode": completed.returncode,
        "stdout": completed.stdout.strip(),
        "stderr": completed.stderr.strip(),
    }


def _git_info() -> dict[str, str]:
    values: dict[str, str] = {}
    for key, command in {
        "commit": ["git", "rev-parse", "HEAD"],
        "branch": ["git", "branch", "--show-current"],
        "status": ["git", "status", "--short"],
    }.items():
        result = _command(command)
        values[key] = str(result.get("stdout", "unknown"))
    return values


def collect_probe() -> dict[str, Any]:
    import cv2
    import numpy
    import scipy
    import sim_swim

    return {
        "pipeline_name": "cs10_runtime_probe",
        "created_at": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(),
        "host": {
            "platform": platform.platform(),
            "python_executable": sys.executable,
            "python_version": sys.version,
            "logical_cpu_count": os.cpu_count(),
            "glibc": platform.libc_ver(),
        },
        "system_commands": {
            "lscpu": _command(["lscpu"]),
            "free_h": _command(["free", "-h"]),
            "ldd_version": _command(["ldd", "--version"]),
            "nvidia_smi": _command(
                [
                    "nvidia-smi",
                    "--query-gpu=name,driver_version,memory.total",
                    "--format=csv,noheader",
                ]
            ),
        },
        "thread_environment": {name: os.environ.get(name) for name in THREAD_VARIABLES},
        "packages": {name: importlib.metadata.version(name) for name in PACKAGES},
        "imports": {
            "sim_swim": str(sim_swim.__file__),
            "cv2": cv2.__version__,
            "numpy": numpy.__version__,
            "scipy": scipy.__version__,
        },
        "numpy_build_config": numpy.__config__.CONFIG,
        "git": _git_info(),
    }


def main(argv: list[str] | None = None) -> Path:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    manifest = collect_probe()
    manifest_path = args.output_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    (args.output_dir / "run.log").write_text(
        "\n".join(
            [
                "pipeline_name=cs10_runtime_probe",
                f"created_at={manifest['created_at']}",
                f"python={manifest['host']['python_executable']}",
                f"manifest={manifest_path}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(manifest_path)
    return manifest_path


if __name__ == "__main__":
    main()
