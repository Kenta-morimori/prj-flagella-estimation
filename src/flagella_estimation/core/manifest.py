from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

from flagella_estimation.core.utils.git import GitInfo
from flagella_estimation.core.utils.time import RunTime


def write_manifest(
    path: str | Path,
    git: GitInfo,
    rt: RunTime,
    input_info: dict[str, Any],
    outputs: dict[str, Any],
) -> None:
    manifest = {
        "git": asdict(git),
        "run_time": asdict(rt),
        "input": input_info,
        "outputs": outputs,
    }
    Path(path).write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
