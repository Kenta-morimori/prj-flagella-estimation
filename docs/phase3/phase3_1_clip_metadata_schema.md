# Phase 3.1 common clip metadata schema

この文書は Issue #127 の実動画・擬似動画共通 clip / metadata schema の正本である。

## Scope

対象は Phase 3 が Phase 4 へ渡す個体 clip と metadata の契約である。

ここでは次を決める。

- 実動画 detection / tracking 経路と Phase 2 擬似動画 GT passthrough 経路の共通出力
- `source_video` / `track` / `clip` / `frame` の単位
- split leakage を防ぐ `group_key`
- label と provenance の責務
- 最小 fixture と library-level contract test

ここでは次を決めない。

- 実動画 detection model の採用
- 1 clip の採用時間長と必要独立 run 数
- 学習 dataset に混ぜてよい物理条件変更
- `n_flagella=4` 以上の training candidate 復帰

## Processing Modes

| mode | input | detection | tracking | label |
| --- | --- | --- | --- | --- |
| `detection_tracking` | 実顕微鏡動画 | required | required | 原則 unavailable |
| `gt_passthrough` | Phase 2 擬似動画 | skipped | GT track を使う | Phase 2 GT `n_flagella` |

どちらの mode でも Phase 3 出力 schema は同じにする。違いは `processing_mode` と `label_source` / provenance で表現する。

## Data Units

| unit | role | independent sample扱い |
| --- | --- | --- |
| `source_video` | 元動画または render 済み擬似動画 | 単独では学習 sample ではない |
| `track` | 1個体の時系列軌跡 | split group の基本単位 |
| `clip` | 固定長 frame sequence | Phase 4 の入力候補 |
| `frame` | clip 内の1時刻 | 独立 sample ではない |

同一 `run_id` / `source_video_id` / `track_id` 由来の `clip` は train / validation / test をまたがせない。Phase 3 metadata では `track.group_key` を必須とし、Phase 4 split はこの値を使う。

## Required IDs

| field | required | responsibility |
| --- | --- | --- |
| `dataset_id` | yes | Phase 3 出力 dataset の論理 ID |
| `source_video.source_video_id` | yes | 元動画を識別する |
| `track.track_id` | yes | source video 内または dataset 内の track を識別する |
| `track.group_key` | yes | split leakage 防止の group |
| `clip.clip_id` | yes | Phase 4 入力候補 clip を識別する |
| `frames[].frame_id` | yes | clip 内 frame を識別する |
| `provenance.run_id` | pseudo only | Phase 2 simulation run を識別する |
| `provenance.render_id` | pseudo render variation when available | 同一 raw run 由来 render を束ねる |

`clip_id` や `render_id` だけを split group に使ってはいけない。

## Metadata Schema

機械可読 schema は `schemas/phase3_clip_metadata.schema.json` に置く。最小例は `examples/phase3/clip_metadata_minimal.json` である。

top-level required fields:

- `schema_version`
- `dataset_id`
- `source_video`
- `processing_mode`
- `provenance`
- `track`
- `clip`
- `normalization`
- `frames`
- `labels`
- `qc`

## Coordinate And Unit Conventions

- pixel 座標は source video の左上原点とする。
- bbox は `[x, y, width, height]` の `bbox_xywh_px` で表す。
- `t_s` / `duration_s` は秒単位。
- `frame_rate_hz` は Hz。
- crop 後画像列のサイズは `normalization.crop_size_px = [width, height]` とする。
- scale が未確定または未使用の場合は `scale_mode=none` とし、`pixel_size_um` は省略してよい。

## Labels

擬似動画では `labels.n_flagella` を Phase 2 GT として保持する。

実動画では `labels.n_flagella` は原則 `null` とし、`label_source=unavailable` にする。手動ラベルが将来入る場合は `label_source=manual` とする。

## QC

Phase 3 は clip ごとに少なくとも次を記録する。

- `qc.status`: `pass`, `review`, `fail`
- `qc.exclusion_reason`: 除外理由。通過時は `null`

実動画では detection confidence / tracking gap、擬似動画では GT bbox / render availability など、mode 固有の詳細は追加 field で拡張してよい。ただし共通 required field は変えない。

## Implementation Boundary

Issue #127 の実装は schema / fixture / contract test までとする。実際の detection、tracking、crop生成 CLI、短時間 clip 評価、learning curve は #6 / #129 に分離する。

