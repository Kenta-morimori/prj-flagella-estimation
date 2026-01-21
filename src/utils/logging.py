from __future__ import annotations

import logging
from pathlib import Path


def setup_logging(log_path: str | Path) -> logging.Logger:
    logger = logging.getLogger("flagella")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()

    fmt = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # File
    fh = logging.FileHandler(Path(log_path), encoding="utf-8")
    fh.setFormatter(fmt)
    fh.setLevel(logging.INFO)

    # Console
    ch = logging.StreamHandler()
    ch.setFormatter(fmt)
    ch.setLevel(logging.INFO)

    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger
