# Phase 3 issue map

この文書は Phase 3 / 4 MVP へ進むための着手順を管理する短い issue map である。

## Fixed Decisions

2026-07-23 のユーザー判断に基づき，次を固定する。

- Phase 3 handoff baseline は Phase 2 dataset v1 の RUN固定 `n_flagella=1,2,3` とする。
- Phase 3 / 4 MVP の標準 clip duration は `0.5 s` とする。`0.25 s` / `1.0 s` は比較条件とする。
- `n_flagella=4` は v1/MVP では diagnostic-only とし，v2 で再検討する。`n_flagella>=5` は対象外とする。
- TUMBLE / Brownian / torque variation / model変更は MVP baseline 外とする。
- RUN-TUMBLE は v2以降の別 dataset とし，論文例の `1200 s` scale ではなく dataset 用短縮 profile を定義する。
- render variation は軽い観測 augmentation のみ training に含める。
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
| 4 | #6 | first implementation | Phase 3 pipeline を実装する | Phase 2 pseudo GT passthrough CLI，`.npy` clip，#127 metadata JSONL，grouped split / QC summaries |

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

#127 は schema と contract test までで完了した。#6 の最初の実装は Phase 2 pseudo GT passthrough に限定し，実動画 detection / tracking は #8 / #9 後に分離する。重い動画生成，長時間 simulation，clip時間長の最終評価，training dataset への条件混在判断は #129 / #128 / Phase 4 へ残す。
