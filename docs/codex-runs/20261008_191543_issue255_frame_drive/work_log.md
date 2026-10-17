# Issue #255：frame/motor修正候補の実装

- 実行target: execution:cs10 / cs10_user_run。ローカルでの実装・targeted test・dry-runまで。simulation・cs10操作なし。
- 作業branch: codex/issue255-frame-drive。PR #261のhead 555cc6aを基点に修正。
- ユーザーのnf3～5各1条件指定に合わせ、両修正Cの計3条件へ縮小。A/Bは設定と単体検証で分離し、動的比較は今回のcampaignに含めない。
- frameは正規化・projection・外積・layer重心の解析的逆伝播による全位置勾配。有限差分勾配、force/torque保存、剛体変換を検証。
- motorは実際のdiffusive segment weightsから計算した力へ、全べん毛ビーズ上のゼロ合力・所定全トルクを満たす最小ノルム補正。縮退は明示エラー。力の増大まで解消したという主張はしない。
- 新規診断はpre-stepの同一座標で成分ごとの総力・総トルクとmotor correction/support/solver/fallbackを記録。online summaryの極値に評価時刻を付加。
- 回帰検証中にbody診断CSVへの余分な列挿入と既存2015 manifestの既定値キー追加を検出し修正。最終状態の全体検証を実施。
- 1tau/1sのCLI dry-runで各3条件、preflightとoutput分離、実効3 workers。outputs/2026-10-08/191543/issue255_frame_drive_preparation/にplan・manifest・run.logを保存。
- hook stiffness・geometry・排除体積・90°拘束・共通QC閾値は変更なし。queue #21は操作しない。既存reservationを新commitへ置換しない。
- 新規ADR 0025はProposed。束化と安定遊泳の改善、単独修正の寄与、canonical採択は未判断。

最終検証: full pytest 984 PASS（200.67 s）、ruff check/format check、git diff --check PASS。local reviewは実装・実行準備についてPASS。
