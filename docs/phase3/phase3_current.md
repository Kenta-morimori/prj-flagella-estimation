# Phase 3 Current

このファイルは，Phase 3 作業の入口として読む短い現在地ドキュメントである。

## Goal

Phase 3 の目的は，実顕微鏡動画と Phase 2 擬似動画を，Phase 4 で共通に使える個体clipとmetadataへ変換することである。

## Current Status

Phase 3 は本格実装前である。Phase 2 dataset v1 の RUN固定 `n_flagella=1,2,3` が training candidate として handoff されている。

現在の主対象:

- #126: 2D投影後にも `n_flagella` 差が残るか確認する。
- #127: 実動画 detection 経路と擬似動画 GT passthrough 経路を，共通clip / metadata schemaへ収束させる。
- #129: 1 clip の時間長と必要な独立run数を決める。
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

1. #127 で共通clip / metadata schemaを文書化する。
2. #129 で clip duration と grouped split 評価を設計する。
3. #6 で最小pipeline実装に入る。
