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

## Phase 3.2: clip duration / dataset mixing / pipeline prep

### P3-2-001: Issues #129 / #128 / #6 Phase 3 schema後続設計を整理する

- status: complete
- source issues:
  - `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/129`
  - `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/128`
  - `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/6`
- branch: `docs/phase3-roadmap-housekeeping-129-128`
- goal: #127 / PR #142 merge後の roadmap housekeeping を行い，clip duration，dataset mixing / versioning，最小 common clip pipeline 実装準備を #127 schema と接続する。
- result:
  - #127 が closed / merged 済みであり，active PR がないことを `docs/ROADMAP.md`，`docs/phase3/phase3_current.md`，`docs/phase3/phase3_issue_map.md` に反映した。
  - #129 の用語，`track.group_key` による grouped split，`0.25 s` / `0.5 s` / `1.0 s` と non-overlap / overlap の評価設計を `docs/phase3/phase3_2_clip_duration_run_count.md` に整理した。
  - #128 の観測augmentation，独立run，domain variation，dataset version変更の分類と，`model_id` / `render_id` / `dataset_version` / `group_key` の扱いを `docs/phase3/phase3_3_dataset_mixing_versioning.md` に整理した。
  - #6 の最小実装を擬似動画 GT passthrough から始める作業分解，入出力，test境界を `docs/phase3/phase3_4_common_clip_pipeline_plan.md` に整理した。
  - 追加判断として，MVP default clip duration を `0.5 s`，torque variation を diagnostic-only，render variation を軽い augmentation のみ，Brownian を当面対象外，RUN-TUMBLE を v2以降の短縮 profile，`n_flagella=4` をv2再検討対象として固定した。
  - 重い learning curve / pilot clip dataset 実行は未実行とし，draft exact command，出力先，判断ポイントを docs に残した。
- acceptance criteria:
  - [x] #127 の closed / merged 状態が docs と整合している。
  - [x] clip duration / source duration / simulation run / track / window の用語が整理されている。
  - [x] grouped split と leakage 防止規則が #127 schema の `group_key` と接続されている。
  - [x] `0.25 s` / `0.5 s` / `1.0 s`，non-overlap / overlap の評価設計がある。
  - [x] augmentation / domain variation / dataset version 規則が #127 schema provenance と接続されている。
  - [x] #6 の最小 GT passthrough pipeline 実装境界が整理されている。
- verification:
  - `git diff --check`
  - `uv run ruff format --check .`
  - `uv run ruff check .`
- docs:
  - `docs/ROADMAP.md`
  - `docs/phase3/phase3_current.md`
  - `docs/phase3/phase3_issue_map.md`
  - `docs/phase3/phase3_1_clip_metadata_schema.md`
  - `docs/phase3/phase3_2_clip_duration_run_count.md`
  - `docs/phase3/phase3_3_dataset_mixing_versioning.md`
  - `docs/phase3/phase3_4_common_clip_pipeline_plan.md`
