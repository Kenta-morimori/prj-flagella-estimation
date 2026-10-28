# Issue263 nf3・1.5秒比較の完了評価

2026-10-28ユーザーはcs10操作と完了確認を明示許可。#39 succeeded、18:19 JST全完走。source commit236d6638を固定し、全5条件375000step/70files/全SHA照合後にrun_summary先読、共通QC5/5PASS、standard/qualitative3D/2D replayと共通bundle診断を実施。

raw archive終端に4.9ps差の保存点があり差分速度を歪めたため、共通bundle解析のみでdt0.1%以下かつ座標1e-6b以下の重複を検証・統合。原archive/QC/シミュレーションは不変。targeted14/full1096PASS。定量表は修正後解析を使い、修正前解析もローカルに保持。

モデル採択/merge/新規予約なし。B/M/C/G/Rの全時間前方変位は負で、C/Rに小さい根元横距離と大きいmotor横トルク。因果未確定/周期不足を明記。Issue/PRは重要項目の既存コメントへ反映する。
