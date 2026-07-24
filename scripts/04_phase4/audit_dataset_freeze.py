#!/usr/bin/env python3
"""Audit a Phase 3 common clip dataset against the Phase 4 freeze policy."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from flagella_estimation.phase4.freeze_workflow import (
    load_freeze_audit_config,
    run_freeze_audit,
)


def _split_config_key(argv: list[str]) -> tuple[Path | None, list[str]]:
    rest: list[str] = []
    config: Path | None = None
    for item in argv:
        if item.startswith("config="):
            config = Path(item.split("=", 1)[1])
        else:
            rest.append(item)
    return config, rest


def main(argv: list[str] | None = None) -> None:
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    config_from_key, parser_argv = _split_config_key(raw_argv)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("overrides", nargs="*", help="Optional KEY=VALUE overrides")
    args = parser.parse_args(parser_argv)
    if config_from_key is not None and args.config is not None:
        parser.error("Use either config=PATH or --config PATH (not both)")
    cfg = load_freeze_audit_config(config_from_key or args.config, args.overrides)
    output_dir, audit = run_freeze_audit(cfg)
    print(output_dir)
    if audit.status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
