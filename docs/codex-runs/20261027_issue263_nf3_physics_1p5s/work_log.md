# nf3・1.5秒の5仮説比較

ユーザーは実行条件を了承しcs10操作・予約を許可、時間を1.5秒と確定した。M/C/Gのopt-in実装、B/R対照、共通preflightの初期力記録、2段階campaign/parallel YAMLを追加。Macでは短いlibrary testとfull regression、静的preview/dry-runのみ。

Targeted84/full1094PASS。Bは0017実装と56配列一致。固定commitを確定後、既存dispatcherを利用し先行1tauを予約する。本比較をゲート前に無条件enqueueしない。
