# Phase 3 Current

このファイルは，Phase 3 作業の入口として読む短い現在地ドキュメントである。

## Goal

Phase 3 の目的は，実顕微鏡動画と Phase 2 擬似動画を，Phase 4 で共通に使える個体clipとmetadataへ変換することである。

## Current Status

Phase 3 は common clip schema の固定を #127 / PR #142 で完了し，clip duration / dataset mixing / 最小 pipeline 実装準備へ進んでいる。Phase 2 dataset v1 の RUN固定 `n_flagella=1,2,3` が training candidate として handoff されている。

現在の主対象:

- #127: closed / merged。実動画 detection 経路と擬似動画 GT passthrough 経路を，共通clip / metadata schemaへ収束させた。schema 正本は `docs/phase3/phase3_1_clip_metadata_schema.md`，機械可読schemaは `schemas/phase3_clip_metadata.schema.json`。
- #129: 1 clip の時間長と必要な独立run数を決める。
- #128: 学習datasetへ混ぜてよい条件変更を整理する。
- #6: 共通clip生成pipelineの実装親Issue。

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

1. #129 の設計に従い，`0.25 s` / `0.5 s` / `1.0 s`，non-overlap / overlap の軽量 fixture と grouped learning curve 実行計画を用意する。
2. #128 の設計に従い，augmentation / domain variation / dataset version 規則を Phase 4 dataset freeze checklist へ接続する。
3. #6 で擬似動画 GT passthrough から最小pipeline実装に入る。実動画 detection / tracking は #8 / #9 の入力条件整理後に本格化する。

## Key references

- Phase 3 issue map: `docs/phase3/phase3_issue_map.md`
- Common clip schema: `docs/phase3/phase3_1_clip_metadata_schema.md`
- Clip duration and run count design: `docs/phase3/phase3_2_clip_duration_run_count.md`
- Dataset mixing and versioning rules: `docs/phase3/phase3_3_dataset_mixing_versioning.md`
- Minimal pipeline implementation plan: `docs/phase3/phase3_4_common_clip_pipeline_plan.md`
- Machine-readable schema: `schemas/phase3_clip_metadata.schema.json`
- Minimal fixture: `examples/phase3/clip_metadata_minimal.json`
