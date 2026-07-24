"""Grouped split helpers for Phase 3 datasets."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Mapping
from dataclasses import dataclass


@dataclass(frozen=True)
class SplitAssignment:
    group_key: str
    split: str


def assign_grouped_splits(
    group_keys: Iterable[str],
    *,
    group_labels: Mapping[str, int] | None = None,
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

    if group_labels is not None:
        missing = [group_key for group_key in groups if group_key not in group_labels]
        if missing:
            raise ValueError(f"group_labels missing group keys: {missing}")
        by_label: dict[int, list[str]] = defaultdict(list)
        for group_key in groups:
            by_label[int(group_labels[group_key])].append(group_key)
        assignments: dict[str, str] = {}
        for label_groups in by_label.values():
            assignments.update(
                _assign_ordered_groups(
                    sorted(label_groups),
                    train_fraction=train_fraction,
                    val_fraction=val_fraction,
                )
            )
        return assignments

    return _assign_ordered_groups(
        groups,
        train_fraction=train_fraction,
        val_fraction=val_fraction,
    )


def _assign_ordered_groups(
    groups: list[str],
    *,
    train_fraction: float,
    val_fraction: float,
) -> dict[str, str]:
    """Assign already ordered groups to splits."""

    n_groups = len(groups)
    train_count = min(n_groups, max(1, int(round(n_groups * train_fraction))))
    remaining_after_train = n_groups - train_count
    val_count = min(remaining_after_train, int(round(n_groups * val_fraction)))
    if n_groups >= 3 and val_count == 0:
        val_count = 1
    if n_groups >= 3 and train_count + val_count == n_groups and train_count > 1:
        train_count -= 1
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
