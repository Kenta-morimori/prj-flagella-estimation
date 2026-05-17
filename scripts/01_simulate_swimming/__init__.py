"""Phase2 simulation utilities package."""

from __future__ import annotations

import importlib.util
from pathlib import Path


def _load_impl():
    impl_path = Path(__file__).resolve().parent / "01_simulate_swimming.py"
    spec = importlib.util.spec_from_file_location("scripts_01_swim_impl", impl_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to load module from: {impl_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_impl = _load_impl()
app = _impl.app
main = _impl.main

__all__ = ["app", "main"]
