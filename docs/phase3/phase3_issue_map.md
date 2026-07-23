# Phase 3 issue map

この文書は Phase 3 / 4 MVP へ進むための着手順を管理する短い issue map である。

## Fixed Decisions

- Phase 3 handoff baseline は Phase 2 dataset v1 の RUN固定 `n_flagella=1,2,3` とする。
- `n_flagella=4` は diagnostic-only とし，`n_flagella>=5` は対象外とする。
- TUMBLE / Brownian / torque variation / model変更は baseline 外とする。
- 同一 `run_id` / `source_video_id` / `track_id` 由来 clip は train / validation / test をまたがせない。

## Map

```text
Phase 3 / 4 MVP
└─ #6 Phase 3 common clip pipeline
   ├─ #127 common clip / metadata schema
   ├─ #129 clip duration and independent run count
   ├─ #128 dataset mixing / augmentation / versioning rules
   └─ #6 implementation PRs
```

## Current Order

| order | issue | state | purpose | deliverable |
| ---: | --- | --- | --- | --- |
| 1 | #127 | closed; PR #142 merged | 実動画・擬似動画の共通 clip / metadata schema を固定する | `docs/phase3/phase3_1_clip_metadata_schema.md`, `schemas/phase3_clip_metadata.schema.json`, 最小 fixture |
| 2 | #129 | in design | clip時間長と必要独立run数を評価する | `docs/phase3/phase3_2_clip_duration_run_count.md`，grouped split 評価計画，learning curve 実行案 |
| 3 | #128 | in design | 学習datasetへ混ぜてよい条件変更を分類する | `docs/phase3/phase3_3_dataset_mixing_versioning.md`，augmentation / domain variation / dataset version規則 |
| 4 | #6 | implementation prep | Phase 3 pipeline を実装する | `docs/phase3/phase3_4_common_clip_pipeline_plan.md`，GT passthrough / crop / metadata 出力 |

## Supporting Issues

| issue | role | blocks MVP baseline? |
| --- | --- | --- |
| #8 | 実顕微鏡動画の取得・入力条件 | 実動画経路の本格実装前に必要 |
| #9 | 実菌体長さ特徴と scale normalization 根拠 | scale normalization 採用時に必要 |
| #17 | 擬似顕微鏡動画の描画条件 | pseudo-video domain tuning 時に必要 |
| #124 | `n>=4` の物理改善 | no |
| #69 | TUMBLE 状態 | no |
| #15 | Brownian 項 | no |

## Implementation Boundary

#127 は schema と contract test までで完了した。重い動画生成，長時間 simulation，clip時間長の採用判断，training dataset への条件混在判断は #129 / #128 / #6 へ残す。
