"""Backward-compatible wrapper for scripts.01_simulate_swimming.plot_motor_scale_collapse_heatmap."""

from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Any


def _load_impl():
    impl_path = (
        Path(__file__).resolve().parent
        / "01_simulate_swimming"
        / "plot_motor_scale_collapse_heatmap.py"
    )
    spec = importlib.util.spec_from_file_location(
        "scripts_plot_heatmap_impl", impl_path
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to load module from: {impl_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_impl = _load_impl()
main = _impl.main


def __getattr__(name: str) -> Any:
    return getattr(_impl, name)


if __name__ == "__main__":
    main()
