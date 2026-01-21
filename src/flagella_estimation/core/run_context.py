from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from flagella_estimation.core.io.paths import OutputPaths, make_output_dirs
from flagella_estimation.core.manifest import write_manifest
from flagella_estimation.core.utils.git import GitInfo, require_committed_state
from flagella_estimation.core.utils.logging import setup_logging
from flagella_estimation.core.utils.time import RunTime, now_jst


@dataclass(frozen=True)
class RunContext:
    git: GitInfo
    run_time: RunTime
    out: OutputPaths
    logger_name: str = "flagella"


def init_run(base_dir: str | Path, input_info: dict[str, Any]) -> RunContext:
    git = require_committed_state()
    rt = now_jst()
    out = make_output_dirs(base_dir=base_dir, date=rt.date, time=rt.time)

    log_path = out.root / "run.log"
    logger = setup_logging(log_path)

    logger.info("Run start")
    logger.info("git.commit=%s", git.commit)
    logger.info("git.branch=%s", git.branch)
    logger.info("output_dir=%s", str(out.root))

    write_manifest(
        path=out.root / "manifest.json",
        git=git,
        rt=rt,
        input_info=input_info,
        outputs={
            "root": str(out.root),
            "tracking_dir": str(out.tracking_dir),
        },
    )

    return RunContext(git=git, run_time=rt, out=out)
