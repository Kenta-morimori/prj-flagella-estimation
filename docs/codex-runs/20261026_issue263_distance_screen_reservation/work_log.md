# Issue #263：距離修正screenの予約と残タスク

- ユーザーはcs10予約を許可し、残タスクと事前判断点の整理を依頼した。
- 予約範囲は準備済みの固定0017fbb・nf6の4us/2us・1tauのみ。1秒stageのQC/replay前予約はしない。
- SSH接続、queue/dispatcher/NAS確認。queue unpaused・他active0。既存dispatcherを継続利用した。
- GitHubから固定commitをfetchし、一時worktreeでsimulationなしremote dry-runと4ファイルのlocal/remote SHA一致を確認した。予備のlocal bundleは作成したが、転送・使用していない。
- 12:36:29 JSTにpriority0で#37をenqueue、12:36:53 JSTに同じ固定commit・2条件のlaunchを確認。job running/failed0、実効2workers・数値thread各1。
- 出力root・予約receipt・remote preflight・launchをローカルへ保存。稼働中のモデル・予約・dispatcherを変更しない。
- runbookに残タスクの順序、1秒比較後の採用dt/26条件/観測周期数の判断、マージ前の位置復元処理削除を整理した。
- Issueは既存の準備要約を予約状況へ更新し、詳細はrunbook/PRを参照させる。新規の定期コメントは追加しない。
- 自動確認なしの既存指示との関係を質問し、ユーザーが今回13:40 JSTの1回確認と完走後のQC/replayを許可した。既存automationをこの1回だけへ更新し、以降の自動再確認を行わない。科学的採択・mergeなし。
