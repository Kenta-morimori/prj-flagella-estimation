"""Opt-in telemetry from the production segment force pass; never a force or gate."""

from __future__ import annotations

import numpy as np

CONTACT_GROUPS = ("bf", "inter_ff", "self_ff")
CONTACT_SUFFIXES = (
    "candidate_pairs",
    "active_pairs",
    "min_distance_um",
    "max_pair_force_N",
    "max_bead_force_N",
    "witness_segment_a",
    "witness_segment_b",
    "witness_s",
    "witness_t",
    "first_active_t_s",
    "active_time_s",
    "first_segment_a",
    "first_segment_b",
)
CONTACT_COLUMNS = ["segment_contact_available", "segment_contact_recorded_steps"] + [
    f"segment_contact_{g}_{suffix}"
    for g in CONTACT_GROUPS
    for suffix in CONTACT_SUFFIXES
]
CONTACT_POST_STEP_COLUMNS = {"segment_contact_recorded_steps"} | {
    f"segment_contact_{g}_active_time_s" for g in CONTACT_GROUPS
}


class SegmentContactStep:
    def __init__(self, groups: np.ndarray, bead_count: int):
        self.groups = groups
        self.candidates = np.bincount(groups, minlength=3)
        self.active = np.zeros(3, dtype=int)
        self.minimum = np.full(3, np.inf)
        self.max_force = np.zeros(3)
        self.witness = np.full((3, 4), -1.0)
        self.forces = np.zeros((3, bead_count, 3))

    def distance(self, index, segment_a, segment_b, s, t, distance):
        group = int(self.groups[index])
        if distance < self.minimum[group]:
            self.minimum[group] = distance
            self.witness[group] = (segment_a, segment_b, s, t)

    def applied_force(self, index, beads, s, t, force):
        group = int(self.groups[index])
        self.active[group] += 1
        self.max_force[group] = max(self.max_force[group], float(np.linalg.norm(force)))
        i, j, k, ell = beads
        self.forces[group, i] += (1.0 - s) * force
        self.forces[group, j] += s * force
        self.forces[group, k] -= (1.0 - t) * force
        self.forces[group, ell] -= t * force


class SegmentContactTracker:
    def __init__(self, springs, pairs, bead_is_body, bead_flag_ids):
        body = np.all(bead_is_body[springs], axis=1)
        owner = np.max(bead_flag_ids[springs], axis=1)
        bf = body[pairs[:, 0]] ^ body[pairs[:, 1]]
        inter = (~bf) & (owner[pairs[:, 0]] != owner[pairs[:, 1]])
        self.groups = np.where(bf, 0, np.where(inter, 1, 2))
        self.bead_count = len(bead_is_body)
        self.first = np.full(3, -1.0)
        self.first_segments = np.full((3, 2), -1.0)
        self.duration = np.zeros(3)
        self.steps = 0

    def new_step(self):
        return SegmentContactStep(self.groups, self.bead_count)

    def metrics(self, step, *, evaluation_t_s, dt_s):
        # Frozen-state evaluation must not advance contact history.
        if dt_s > 0:
            active = step.active > 0
            new = (self.first < 0) & active
            self.first[new] = evaluation_t_s
            self.first_segments[new] = step.witness[new, :2]
            self.duration[active] += dt_s
            self.steps += 1
        result = {
            "segment_contact_available": 1,
            "segment_contact_recorded_steps": self.steps,
        }
        for i, group in enumerate(CONTACT_GROUPS):
            values = (
                int(step.candidates[i]),
                int(step.active[i]),
                float(step.minimum[i] * 1e6) if step.candidates[i] else -1.0,
                float(step.max_force[i]),
                float(np.linalg.norm(step.forces[i], axis=1).max(initial=0)),
                *step.witness[i],
                float(self.first[i]),
                float(self.duration[i]),
                *self.first_segments[i],
            )
            result.update(
                {
                    f"segment_contact_{group}_{k}": v
                    for k, v in zip(CONTACT_SUFFIXES, values, strict=True)
                }
            )
        return result
