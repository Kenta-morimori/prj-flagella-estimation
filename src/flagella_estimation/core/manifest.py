from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

from flagella_estimation.core.utils.git import GitInfo
from flagella_estimation.core.utils.time import RunTime


@dataclass(frozen=True)
class Manifest:
    git: dict[str, Any]
    run_time: dict[str, Any]
    input: dict[str, Any]
    outputs: dict[str, Any]


def write_manifest(path: str | Path, git: GitInfo, rt: RunTime, input_info: dict, outputs: dict) -> None:
    m = Manifest(
        git={"commit": git.commit, "commit_short": git.commit_short, "branch": git.branch},
        run_time={"date": rt.date, "time": rt.time},
        input=input_info,
        outputs=outputs,
    )
    Path(path).write_text(json.dumps(asdict(m), ensure_ascii=False, indent=2), encoding="utf-8")
