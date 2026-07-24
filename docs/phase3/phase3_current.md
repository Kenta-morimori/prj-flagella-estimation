# Phase 3 Current

このファイルは，Phase 3 作業の入口として読む短い現在地ドキュメントである。

## Goal

Phase 3 の目的は，実顕微鏡動画と Phase 2 擬似動画を，Phase 4 で共通に使える個体clipとmetadataへ変換することである。

## Current Status

Phase 3 は common clip schemaを#127 / PR #142，pseudo GT passthroughを#6 / PR #144で実装済みである。Phase 2 dataset v1 の RUN固定 `n_flagella=1,2,3` がtraining candidateとしてPhase 4 loader / baseline / learning curve / freeze gateまで接続されている。

現在の主対象:

- #127: closed / merged。実動画 detection 経路と擬似動画 GT passthrough 経路を，共通clip / metadata schemaへ収束させた。schema 正本は `docs/phase3/phase3_1_clip_metadata_schema.md`，機械可読schemaは `schemas/phase3_clip_metadata.schema.json`。
- #129: `0.5 s` defaultとgrouped learning curveは完了し，pseudo-v1で`k=4`を採用する範囲を判断する。
- #128: closed / PR #152 merged。Phase 4 machine-readable freeze gateへ接続した。
- #6: 共通clip生成pipelineの実装親Issue。Phase 2 擬似動画 GT passthrough は実装済みで，実動画 detection / tracking は #8 / #9 後に進める。

MVP 固定方針（2026-07-23 のユーザー判断に基づく）:

- Phase 3 / 4 MVP の標準 clip duration は `0.5 s` とする。`0.25 s` / `1.0 s` は比較条件に残す。
- torque variation は MVP training baseline に混ぜず，diagnostic / robustness dataset として分離する。
- Brownian は当面含めない。
- RUN-TUMBLE は v2 以降で，論文の状態遷移構造を使いつつ短縮 profile を別定義する。
- render variation は軽い観測 augmentation のみ training に含める。
- `n_flagella=4` は v1/MVP では diagnostic-only のまま，v2 で再検討する。

## Input Paths

Phase 2擬似動画:

- raw simulation: `outputs/phase2_multi_run/flagella_count_behavior_v1`
- analysis dataset: `outputs/phase2_analysis/flagella_count_behavior/datasets/v1`
- current training candidate: `n_flagella=1,2,3`

実顕微鏡動画:

- 入力条件は #8 で整理する。
- 菌体長さ特徴は #9 で整理する。

## Design Boundary

Phase 3 の出力schemaは，detectionあり・GT passthroughのどちらでも同じにする。

最低限保持するID:

- `source_video_id`
- `track_id`
- `clip_id`
- `frame_id`
- `run_id` または simulation provenance
- `group_key`

dataset splitでは，同一 `run_id` / `source_video_id` / `track_id` 由来のclipを train / validation / test にまたがせない。

## Next Actions

1. #129で`k=4`をpseudo-v1 MVP lower boundとして採用するか判断する。
2. #128 freeze gateを今後のtraining / learning curve前にも継続実行する。
3. #145 RUN-TUMBLE dataset v2はProject `TODO`のまま後回しにする。
4. 実動画 detection / tracking は #8 / #9 の入力条件整理後に本格化する。

## Key references

- Phase 3 issue map: `docs/phase3/phase3_issue_map.md`
- Common clip schema: `docs/phase3/phase3_1_clip_metadata_schema.md`
- Clip duration and run count design: `docs/phase3/phase3_2_clip_duration_run_count.md`
- Dataset mixing and versioning rules: `docs/phase3/phase3_3_dataset_mixing_versioning.md`
- Minimal pipeline implementation plan: `docs/phase3/phase3_4_common_clip_pipeline_plan.md`
- Machine-readable schema: `schemas/phase3_clip_metadata.schema.json`
- Minimal fixture: `examples/phase3/clip_metadata_minimal.json`
