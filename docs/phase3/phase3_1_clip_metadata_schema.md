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

機械可読 schema では，`processing_mode=gt_passthrough` または `source_video.source_kind` が pseudo 系の場合に，non-null の `provenance.run_id` を要求する。

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
- `source_video.frame_count`，`source_video.codec_fourcc`，`source_video.file_size_bytes` は実動画入力の再現性確認用の任意 field である。
- crop 後画像列のサイズは `normalization.crop_size_px = [width, height]` とする。
- scale が未確定または未使用の場合は `scale_mode=none` とし、`pixel_size_um` は省略してよい。
- 実動画 detection では，bbox だけでなく `frames[].body_axis_angle_rad`，`frames[].body_length_px`，`frames[].body_width_px`，`frames[].detection_confidence` を任意で保持してよい。低コントラストの楕円/桿状候補を後で監査し，`body_axis_align` や `body_length` normalization の妥当性を確認するためである。

## Labels

擬似動画では `labels.n_flagella` を Phase 2 GT として保持する。

実動画では `labels.n_flagella` は原則 `null` とし、`label_source=unavailable` にする。手動ラベルが将来入る場合は `label_source=manual` とする。

機械可読 schema では，`label_source=unavailable` のとき `n_flagella=null`，`label_source=phase2_gt` または `manual` のとき `n_flagella` は0以上の整数であることを要求する。

## QC

Phase 3 は clip ごとに少なくとも次を記録する。

- `qc.status`: `pass`, `review`, `fail`
- `qc.exclusion_reason`: 除外理由。通過時は `null`

実動画では `qc.detection_confidence_min`，`qc.tracking_gap_count`，`qc.notes` を任意で保持してよい。擬似動画では GT bbox / render availability などを同じ `qc.status` / `qc.exclusion_reason` へ畳み込む。共通 required field は変えない。

## Real AVI Initial Analysis Notes

2026-07-23 に `data/20250716data1.avi` と `data/20250716data2.avi` を使って，Otsu + connected components の初期確認を行った。

- 2本とも 512x512, 128 frames で，OpenCV/file 上の reported fps は `20250716data1.avi` が 25.0 Hz，`20250716data2.avi` が 20.0 Hz だった。
- `candidate_count_mean` は約 43-46/frame で，単純な dark connected component は小さい黒点，背景ノイズ，リング状アーチファクトを多く拾った。
- 菌体らしい低コントラスト楕円/桿状構造を後で選別できるように，frame-level の body axis / 長短径 / confidence と clip-level QC を任意 field として残す。
- AVI 本体は `data/` 配下で git 管理外とし，metadata には入力パスや動画要約だけを記録する。

## Issue #127 Close Checklist

- [x] 実動画・擬似動画の共通 required field を固定する。
- [x] Phase 2 GT passthrough と実動画 `label_source=unavailable` の両方を表現できる。
- [x] GT passthrough の `provenance.run_id` と，`label_source` / `n_flagella` の整合性を schema で強制する。
- [x] split leakage 防止用の `track.group_key` を必須にする。
- [x] 実AVI初期分析から，実動画source metadataと検出候補QCに必要な任意fieldを確認する。
- [x] 実 detector / tracker / crop CLI，clip window評価，dataset mixing判断は #6 / #129 / #128 へ分離する。
- [ ] PR上で CI / review gate が通ることを確認する。

## Implementation Boundary

Issue #127 の実装は schema / fixture / contract test までとする。実際の detection、tracking、crop生成 CLI、短時間 clip 評価、learning curve は #6 / #129 に分離する。
