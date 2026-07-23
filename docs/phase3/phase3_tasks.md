# Phase 3 Tasks

このファイルは Phase 3 の採択済みタスクを管理する。

チェックボックスは review PASS 後にのみ更新する。

## Phase 3.1: common clip / metadata schema

### P3-1-001: Issue #127 実動画と擬似動画の共通clip schemaを決める

- status: complete
- source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/127`
- branch: `feature/phase3-127-clip-schema`
- goal: 実動画 detection / tracking 経路と Phase 2 擬似動画 GT passthrough 経路を，Phase 4 が共通利用できる個体clip / metadata schemaへ収束させる。
- result:
  - 共通 schema 正本を `docs/phase3/phase3_1_clip_metadata_schema.md` に追加した。
  - 機械可読 schema を `schemas/phase3_clip_metadata.schema.json` に追加した。
  - 最小 fixture を `examples/phase3/clip_metadata_minimal.json` に追加した。
  - schema required field，`group_key`，Phase 2 GT label，実動画 label unavailable 状態を pytest で固定した。
  - 実装範囲は schema / fixture / contract test までとし，detection / tracking / crop CLI は #6，clip時間長と独立run数は #129，条件混在規則は #128 へ残した。
- acceptance criteria:
  - [x] 実動画・擬似動画の入力差分と共通処理が表で文書化されている。
  - [x] 個体clipとmetadataの必須項目・型・単位が定義されている。
  - [x] detectionを省略してもPhase 3出力schemaが変わらない。
  - [x] dataset split時のgroup keyが保持される。
  - [x] 実装Issueへ分割できる粒度の仕様になっている。
- verification:
  - `uv run pytest tests/test_phase3_clip_metadata_schema.py -q`
- docs:
  - `docs/phase3/phase3_1_clip_metadata_schema.md`
  - `docs/phase3/phase3_issue_map.md`

