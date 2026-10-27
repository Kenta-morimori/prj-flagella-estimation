# nf3・1.5秒の5仮説比較

ユーザーは実行条件を了承しcs10操作・予約を許可、時間を1.5秒と確定した。M/C/Gのopt-in実装、B/R対照、共通preflightの初期力記録、2段階campaign/parallel YAMLを追加。Macでは短いlibrary testとfull regression、静的preview/dry-runのみ。

Targeted84/full1094PASS。Bは0017実装と56配列一致。固定commitを確定後、既存dispatcherを利用し先行1tauを予約する。本比較をゲート前に無条件enqueueしない。

## 短期結果と本比較予約

固定source 236d6638。#38は18:28 JSTに全条件完走、70 files/1461804 bytes全SHA一致、run_summary先読、共通strictQC5/5PASS、固定camera3D/2D各41frame/25fpsを検証・レビュー。日本語fontを描画環境だけで修正し、QC/入力/物理条件は不変。短期の完全周期0でモデル採択はしない。

ユーザーの既存許可範囲で#39を18:34:23 JSTにenqueue、18:34:33 JST正常開始。5条件各375000step、同一固定commit、cs10_qualified5workers、thread各1。他active予約0、既存dispatcher利用、取消/置換/優先度変更なし。短期速度外挿ETA 2026-10-28T17:08:29.747147+09:00。自動監視は追加していない。
