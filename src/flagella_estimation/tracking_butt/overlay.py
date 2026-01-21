from __future__ import annotations

import cv2

from flagella_estimation.tracking_butt.types import ButtEstimate, Detection


class OverlayRenderer:
    def __init__(self, flagella_length: int = 30) -> None:
        self.flagella_length = flagella_length

    def draw(self, frame, detection: Detection, track_id: int, butt: ButtEstimate) -> None:
        """Draw tracking overlays (ellipse/bbox, ids, butt, pseudo flagella)."""
        center_pt = (int(round(detection.cx)), int(round(detection.cy)))
        butt_pt = (int(round(butt.point[0])), int(round(butt.point[1])))

        if detection.is_valid and detection.angle_deg is not None and detection.major and detection.minor:
            ellipse = (
                (float(detection.cx), float(detection.cy)),
                (float(detection.major), float(detection.minor)),
                float(detection.angle_deg),
            )
            cv2.ellipse(
                frame, ellipse, (0, 200, 0), 2
            )
        else:
            x, y, w, h = detection.bbox
            cv2.rectangle(frame, (x, y), (x + w, y + h), (0, 150, 0), 1)

        cv2.circle(frame, center_pt, 3, (0, 255, 255), -1)
        cv2.circle(frame, butt_pt, 4, (0, 0, 255), -1)

        end_pt = (
            int(round(butt.point[0] + butt.flagella_dir[0] * self.flagella_length)),
            int(round(butt.point[1] + butt.flagella_dir[1] * self.flagella_length)),
        )
        cv2.line(frame, butt_pt, end_pt, (255, 0, 255), 2)

        cv2.putText(
            frame,
            f"{track_id}",
            (center_pt[0] + 5, center_pt[1] - 5),
            cv2.FONT_HERSHEY_SIMPLEX,
            0.5,
            (255, 255, 255),
            1,
            cv2.LINE_AA,
        )
