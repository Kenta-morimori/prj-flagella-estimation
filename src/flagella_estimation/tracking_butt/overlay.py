from __future__ import annotations

import cv2

from flagella_estimation.tracking_butt.types import ButtEstimate, Detection


class OverlayRenderer:
    def __init__(
        self,
        flagella_length: int = 30,
        scale_bar_px: float | None = None,
        scale_bar_um: float | None = None,
    ) -> None:
        """擬似べん毛の長さとスケールバー設定を初期化する。

        Args:
            flagella_length: 擬似べん毛を描画する長さ[px]。
            scale_bar_px: スケールバーの長さ[px]。None または 0 以下なら描画しない。
            scale_bar_um: スケールバーに表示する長さ[μm]。None なら px 表記のみ。
        """
        self.flagella_length = flagella_length
        self.scale_bar_px = scale_bar_px if scale_bar_px and scale_bar_px > 0 else None
        self.scale_bar_um = scale_bar_um

    def draw(
        self, frame, detection: Detection, track_id: int, butt: ButtEstimate
    ) -> None:
        """楕円/バウンディングボックス、ID、お尻点、擬似べん毛、スケールバーを描画する。

        Args:
            frame: 描画対象の画像 (BGR)。
            detection: 検知結果。
            track_id: トラックID。
            butt: お尻推定結果。

        Returns:
            None
        """
        center_pt = (int(round(detection.cx)), int(round(detection.cy)))
        butt_pt = (int(round(butt.point[0])), int(round(butt.point[1])))

        if (
            detection.is_valid
            and detection.angle_deg is not None
            and detection.major
            and detection.minor
        ):
            ellipse = (
                (float(detection.cx), float(detection.cy)),
                (float(detection.major), float(detection.minor)),
                float(detection.angle_deg),
            )
            cv2.ellipse(frame, ellipse, (0, 200, 0), 2)
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

        self._draw_scale_bar(frame)

    def _draw_scale_bar(self, frame) -> None:
        """スケールバーを右下に描画する。

        Args:
            frame: 描画対象の画像 (BGR)。

        Returns:
            None
        """
        if self.scale_bar_px is None:
            return

        height, width = frame.shape[:2]
        margin = 20
        length_px = int(round(self.scale_bar_px))
        length_px = max(1, min(length_px, width - 2 * margin))

        end_x = width - margin
        start_x = end_x - length_px
        y = height - margin

        cv2.line(frame, (start_x, y), (end_x, y), (255, 255, 255), 2)
        label = (
            f"{self.scale_bar_um:.2f} um"
            if self.scale_bar_um is not None
            else f"{length_px} px"
        )
        cv2.putText(
            frame,
            label,
            (start_x, y - 8),
            cv2.FONT_HERSHEY_SIMPLEX,
            0.5,
            (255, 255, 255),
            1,
            cv2.LINE_AA,
        )
