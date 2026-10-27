# Issue #263 距離修正1τ screenの評価

- ユーザーのVPN復旧連絡・完了確認/QC/replay用cs10操作許可に基づく。
- cs10はread-only操作とportable読出しのみ。予約・dispatcher・条件・コードの変更なし。
- #37は10/26 13:31:47 JSTに完走済み。固定commit/入力SHA/2条件/stepを照合。
- 31 portableファイルを同期、全SHA一致、run_summary先読、共通QC2/2PASS。
- 共通evaluator/replay/bundle_rotationは実行commitのsource snapshotから実行。
- 条件別/比較replay6本を生成。日本語fontの不足を描画設定で解消し、frame/fps/coverage/SHAと目視を再確認。
- 初期0.04sでは反発力0、旧4usと修正4usは実質同じ、修正4us/2usも近い。回転0.167回、初期過渡の変位は両方負。
- 減速の0.4〜0.8sを含む25tauはユーザーレビュー後。新規予約/実行、torsion比較、モデル採択、mergeなし。
- Issueは既存の予約要約コメントを重要結果のみに編集。詳細をPR/runbook/resultへ保存する。
