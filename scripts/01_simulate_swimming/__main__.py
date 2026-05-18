"""Backward-compatible module entrypoint for `python -m scripts.01_simulate_swimming`."""

from __future__ import annotations

import importlib.util
from pathlib import Path


def _run() -> None:
    impl_path = Path(__file__).resolve().parent / "01_simulate_swimming.py"
    spec = importlib.util.spec_from_file_location("scripts_01_swim_impl", impl_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to load module from: {impl_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.app()


if __name__ == "__main__":
    _run()
