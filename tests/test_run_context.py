from __future__ import annotations

import json
from pathlib import Path

from sim_swim.core.run_context import init_run


def test_init_run_creates_log_and_manifest(tmp_path: Path) -> None:
    ctx = init_run(base_dir=tmp_path, input_info={"config": "conf/sim_swim.yaml"})

    assert ctx.out.root.is_dir()
    assert ctx.out.sim_dir.is_dir()
    assert ctx.out.render_dir.is_dir()
    assert ctx.out.render2d_dir.is_dir()

    log_path = ctx.out.root / "run.log"
    manifest_path = ctx.out.root / "manifest.json"
    assert log_path.is_file()
    assert manifest_path.is_file()

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    assert manifest["input"]["config"] == "conf/sim_swim.yaml"
    assert Path(manifest["outputs"]["sim_dir"]).name == "sim"
