#!/usr/bin/env python3
"""Deprecated alias for the generic Phase 2 shape-stability grid sweep."""

from sim_swim.analysis.sweeps import shape_stability_grid as _impl
from sim_swim.analysis.sweeps.shape_stability_grid import *  # noqa: F401,F403


def __getattr__(name: str):
    return getattr(_impl, name)


def main(argv: list[str] | None = None) -> None:
    _impl.main(argv)


if __name__ == "__main__":
    main()
