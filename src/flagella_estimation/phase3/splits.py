"""Grouped split helpers for Phase 3 datasets."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass


@dataclass(frozen=True)
class SplitAssignment:
    group_key: str
    split: str


def assign_grouped_splits(
    group_keys: Iterable[str],
    *,
    train_fraction: float = 0.67,
    val_fraction: float = 0.22,
) -> dict[str, str]:
    """Assign each group key to exactly one split."""

    groups = sorted(set(group_keys))
    if not groups:
        return {}
    if train_fraction < 0.0 or val_fraction < 0.0:
        raise ValueError("split fractions must be non-negative")
    if train_fraction + val_fraction > 1.0:
        raise ValueError("train_fraction + val_fraction must be <= 1")

    n_groups = len(groups)
    train_count = min(n_groups, int(round(n_groups * train_fraction)))
    val_count = min(n_groups - train_count, int(round(n_groups * val_fraction)))
    assignments: dict[str, str] = {}
    for index, group_key in enumerate(groups):
        if index < train_count:
            split = "train"
        elif index < train_count + val_count:
            split = "val"
        else:
            split = "test"
        assignments[group_key] = split
    return assignments


def assert_no_group_leakage(rows: Iterable[dict[str, str]]) -> None:
    """Raise if the same group key appears in multiple splits."""

    seen: dict[str, str] = {}
    for row in rows:
        group_key = str(row["group_key"])
        split = str(row["split"])
        previous = seen.setdefault(group_key, split)
        if previous != split:
            raise ValueError(
                f"group_key leakage: {group_key} appears in {previous} and {split}"
            )
