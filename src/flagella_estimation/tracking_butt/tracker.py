from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np

from flagella_estimation.tracking_butt.types import Detection, TrackUpdate


@dataclass
class TrackState:
    track_id: int
    last_detection: Detection
    last_frame: int
    last_velocity: Tuple[float, float]
    prev_velocity: Tuple[float, float] | None
    missing_frames: int
    observed_frames: int
    first_frame: int

    @property
    def span_frames(self) -> int:
        return self.observed_frames + self.missing_frames

    def missing_rate(self) -> float:
        denom = self.span_frames
        return self.missing_frames / denom if denom > 0 else 0.0


class Tracker:
    def __init__(self, max_link_distance: float, max_inactive: int = 10) -> None:
        self.max_link_distance = max_link_distance
        self.max_inactive = max_inactive
        self.tracks: Dict[int, TrackState] = {}
        self._next_id = 0

    def _active_tracks(self, frame_idx: int) -> Dict[int, TrackState]:
        return {
            tid: st
            for tid, st in self.tracks.items()
            if frame_idx - st.last_frame <= self.max_inactive
        }

    def step(self, frame_idx: int, detections: List[Detection]) -> List[TrackUpdate]:
        active_tracks = self._active_tracks(frame_idx)

        pairs: list[tuple[float, int, int]] = []
        for d_idx, det in enumerate(detections):
            for t_id, state in active_tracks.items():
                dist = np.hypot(det.cx - state.last_detection.cx, det.cy - state.last_detection.cy)
                if dist <= self.max_link_distance:
                    pairs.append((float(dist), t_id, d_idx))

        pairs.sort(key=lambda x: x[0])
        assigned_tracks: set[int] = set()
        assigned_detections: set[int] = set()
        assignments: list[tuple[int, int]] = []
        for dist, t_id, d_idx in pairs:
            if t_id in assigned_tracks or d_idx in assigned_detections:
                continue
            assignments.append((t_id, d_idx))
            assigned_tracks.add(t_id)
            assigned_detections.add(d_idx)

        for d_idx, det in enumerate(detections):
            if d_idx in assigned_detections:
                continue
            t_id = self._next_id
            self._next_id += 1
            assignments.append((t_id, d_idx))

        updates: list[TrackUpdate] = []
        for t_id, d_idx in assignments:
            det = detections[d_idx]
            state = self.tracks.get(t_id)
            if state is None:
                state = TrackState(
                    track_id=t_id,
                    last_detection=det,
                    last_frame=frame_idx,
                    last_velocity=(0.0, 0.0),
                    prev_velocity=None,
                    missing_frames=0,
                    observed_frames=1,
                    first_frame=frame_idx,
                )
                self.tracks[t_id] = state
                updates.append(
                    TrackUpdate(
                        frame_idx=frame_idx,
                        track_id=t_id,
                        detection=det,
                        vx=0.0,
                        vy=0.0,
                        prev_velocity=None,
                    )
                )
                continue

            frame_gap = max(0, frame_idx - state.last_frame - 1)
            prev_velocity = state.last_velocity
            dt = max(1, frame_idx - state.last_frame)
            vx = (det.cx - state.last_detection.cx) / dt
            vy = (det.cy - state.last_detection.cy) / dt

            state.missing_frames += frame_gap
            state.observed_frames += 1
            state.prev_velocity = prev_velocity
            state.last_velocity = (float(vx), float(vy))
            state.last_detection = det
            state.last_frame = frame_idx

            updates.append(
                TrackUpdate(
                    frame_idx=frame_idx,
                    track_id=t_id,
                    detection=det,
                    vx=float(vx),
                    vy=float(vy),
                    prev_velocity=prev_velocity,
                )
            )
        return updates

    def qc_summary(self) -> Dict[int, Dict[str, float]]:
        summary: Dict[int, Dict[str, float]] = {}
        for t_id, state in self.tracks.items():
            summary[t_id] = {
                "frames": float(state.observed_frames),
                "missing_rate": state.missing_rate(),
            }
        return summary
