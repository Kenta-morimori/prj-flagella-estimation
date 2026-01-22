from __future__ import annotations

import cv2

from flagella_estimation.tracking_butt.types import ButtEstimate, Detection


class OverlayRenderer:
    def __init__(
        self,
        flagella_length: int = 30,
        scale_bar_px: float | None = None,
        scale_bar_um: float | None = None,
        draw_history: bool = True,
        history_length: int = 30,
        hide_history_after: int = 60,
    ) -> None:
        """擬似べん毛の長さとスケールバー設定を初期化する。

        Args:
            flagella_length: 擬似べん毛を描画する長さ[px]。
            scale_bar_px: スケールバーの長さ[px]。None または 0 以下なら描画しない。
            scale_bar_um: スケールバーに表示する長さ[μm]。None なら px 表記のみ。
            draw_history: トラックの軌跡を描画するか。
            history_length: 軌跡として保持する最大フレーム長。
            hide_history_after: このフレーム数より古い履歴は描かない。
        """
        self.flagella_length = flagella_length
        self.scale_bar_px = scale_bar_px if scale_bar_px and scale_bar_px > 0 else None
        self.scale_bar_um = scale_bar_um
        self.draw_history = draw_history
        self.history_length = max(1, history_length)
        self.hide_history_after = max(1, hide_history_after)
        self._history: dict[int, list[tuple[int, float, float]]] = {}

    def draw(
        self,
        frame,
        detection: Detection,
        track_id: int,
        butt: ButtEstimate,
        frame_idx: int | None = None,
    ) -> None:
        """楕円/バウンディングボックス、ID、お尻点、擬似べん毛、スケールバーを描画する。

        Args:
            frame: 描画対象の画像 (BGR)。
            detection: 検知結果。
            track_id: トラックID。
            butt: お尻推定結果。
            frame_idx: 現在のフレーム番号（軌跡描画に使用）。

        Returns:
            None
        """
        center_pt = (int(round(detection.cx)), int(round(detection.cy)))
        butt_pt = (int(round(butt.point[0])), int(round(butt.point[1])))

        if self.draw_history and frame_idx is not None:
            history = self._history.setdefault(track_id, [])
            history.append((frame_idx, detection.cx, detection.cy))
            # 古い履歴を削除
            cutoff = frame_idx - self.hide_history_after
            self._history[track_id] = [
                h for h in history if h[0] >= cutoff
            ][-self.history_length :]
            self._draw_history(frame, self._history[track_id])

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

    def draw_scale_bar(self, frame) -> None:
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

    def _draw_history(
        self, frame, history: list[tuple[int, float, float]]
    ) -> None:
        """トラックの重心履歴を線で描画する。

        Args:
            frame: 描画対象の画像。
            history: (frame_idx, x, y) の履歴リスト（古い順）。
        """
        if len(history) < 2:
            return
        pts = [
            (int(round(x)), int(round(y)))
            for _, x, y in sorted(history, key=lambda t: t[0])
        ]
        for i in range(1, len(pts)):
            cv2.line(frame, pts[i - 1], pts[i], (0, 200, 200), 1)
