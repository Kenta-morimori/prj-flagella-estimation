# ADR 0021: Pending Phase 2モデルの共通評価契約

## Status

Accepted

## Context

新しいPhase 2 model candidateごとに、1τの数値screen、QC、heatmap、replayを個別実装すると、条件軸、出力、gate、長時間評価へのhandoffが揃わない。短時間の数値・物理QCと、長時間archiveを用いる遊泳特徴量評価も目的が異なる。

## Decision

`implementation_status: pending` のPhase 2 profileは、`development_evaluation`契約を持つ。Codexはrepository-discoverable `model-development-evaluation` skillを用い、共通集約器でshort screenまたはlong_durationのcondition coverage、profile、provenance、PASS/FAIL QC、統合summary、window QC、artifact SHA-256、固定camera replayを作成する。

hook angleは診断として保存するが、shape PASS/FAILには用いない。finite、body、hook length、flag、motor diagnosticsは引き続きgateとする。

長時間campaignはshort screenのレビュー後に同じ契約の`long_duration` stageとして定義する。遊泳特徴量評価は、PASSした長時間archiveを入力にする別Issueの責務とする。

`2010_hex_project`のattachment比較では、compound `attachment_pattern` axisに明示slot集合を保存する。C6回転のみを同一視し反射は区別する。同期済みcompleted archiveから、pattern列・べん毛数行のsparse heatmapと固定camera 3D/2D replayを共通集約器で生成する。screen FAIL後に明示許可された連続mainを実行しても、そのmain archiveはdiagnostic-onlyであり、採択根拠にしない。

## Consequences

- Issue専用のmodel-evaluation解析器を新設しない。
- 既存のpending 2010 hex / 2015 project / 2015 paper profileも契約を宣言する。
- 1τ結果だけで長時間束化や特徴量の採否を主張しない。
- long_durationは全conditionのstrict PASS、required artifact、同期後のSHA-256検証を満たすまで後続特徴量評価へ渡さない。
