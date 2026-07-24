#!/usr/bin/env python3
"""Train the Phase 4 common-clip baseline classifier."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

from flagella_estimation.phase4.training import (
    load_baseline_config,
    train_baseline_classifier,
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
    cfg = load_baseline_config(config_from_key or args.config, args.overrides)
    output_dir = train_baseline_classifier(cfg)
    print(output_dir)


if __name__ == "__main__":
    main()
