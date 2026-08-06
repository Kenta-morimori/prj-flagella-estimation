# Phase 3 Common Clip Metadata Contract

この文書は，Phase 3がPhase 4へ渡す個体clipとmetadataのhuman-readable contractである．

machine-readable fieldとvalidation規則の正本は，`schemas/phase3_clip_metadata.schema.json`と対応するtestである．

## Scope

対象は，実顕微鏡動画とPhase 2擬似動画から生成する個体clipとmetadataである．

この文書では次を定める．

* 実動画のdetection / tracking経路と擬似動画のGT passthrough経路の共通出力
* `source_video`，`track`，`clip`，`frame`の単位
* split leakageを防ぐ`group_key`
* labelとprovenanceの責務
* 座標，単位，QCの基本契約

clip時間長，必要独立run数，dataset mixing policyは`phase3_tasks.md`とADR 0015を参照する．

## Processing Modes

| mode                 | input       | detection | tracking    | label                   |
| -------------------- | ----------- | --------- | ----------- | ----------------------- |
| `detection_tracking` | 実顕微鏡動画      | required  | required    | 原則unavailable           |
| `gt_passthrough`     | Phase 2擬似動画 | skipped   | GT trackを使用 | Phase 2 GT `n_flagella` |

どちらのmodeでも出力schemaは同じとする．経路差は`processing_mode`，`label_source`，provenanceで表現する．

## Data Units

| unit           | role               | independent sample扱い |
| -------------- | ------------------ | -------------------- |
| `source_video` | 元動画またはrender済み擬似動画 | 単独では学習sampleではない     |
| `track`        | 1個体の時系列軌跡          | split groupの基本単位     |
| `clip`         | 固定長frame sequence  | Phase 4の入力候補         |
| `frame`        | clip内の1時刻          | 独立sampleではない         |

同一`run_id`，`source_video_id`，`track_id`由来のclipはtrain / validation / testをまたがせない．Phase 4 splitは`track.group_key`を使用する．

## Required IDs

| field                          | required                | responsibility              |
| ------------------------------ | ----------------------- | --------------------------- |
| `dataset_id`                   | yes                     | Phase 3出力datasetの論理ID       |
| `source_video.source_video_id` | yes                     | 元動画を識別する                    |
| `track.track_id`               | yes                     | trackを識別する                  |
| `track.group_key`              | yes                     | split leakageを防ぐ            |
| `clip.clip_id`                 | yes                     | Phase 4入力候補clipを識別する        |
| `frames[].frame_id`            | yes                     | clip内frameを識別する             |
| `provenance.run_id`            | pseudo only             | Phase 2 simulation runを識別する |
| `provenance.render_id`         | available when rendered | 同一raw run由来の描画条件を識別する       |

`clip_id`や`render_id`だけをsplit groupに使ってはいけない．

`processing_mode=gt_passthrough`またはpseudo sourceでは，non-nullの`provenance.run_id`を要求する．

## Metadata Contract

top-level required fieldsは次のとおりである．

* `schema_version`
* `dataset_id`
* `source_video`
* `processing_mode`
* `provenance`
* `track`
* `clip`
* `normalization`
* `frames`
* `labels`
* `qc`

機械可読schemaは`schemas/phase3_clip_metadata.schema.json`，最小例は`examples/phase3/clip_metadata_minimal.json`に置く．

## Coordinate And Unit Conventions

* pixel座標はsource videoの左上原点とする．
* bboxは`[x, y, width, height]`の`bbox_xywh_px`で表す．
* `t_s`と`duration_s`は秒単位とする．
* `frame_rate_hz`はHzとする．
* crop後画像列のサイズは`normalization.crop_size_px = [width, height]`とする．
* scaleが未確定または未使用の場合は`scale_mode=none`とし，`pixel_size_um`は省略できる．
* `source_video.frame_count`，`source_video.codec_fourcc`，`source_video.file_size_bytes`は入力再現性確認用の任意fieldとする．
* 実動画では`frames[].body_axis_angle_rad`，`frames[].body_length_px`，`frames[].body_width_px`，`frames[].detection_confidence`を任意で保持できる．

## Labels

擬似動画では`labels.n_flagella`をPhase 2 GTとして保持し，`label_source=phase2_gt`とする．

実動画では`labels.n_flagella`を原則`null`とし，`label_source=unavailable`とする．手動labelを追加する場合は`label_source=manual`とする．

machine-readable schemaでは次を要求する．

* `label_source=unavailable`: `n_flagella=null`
* `label_source=phase2_gt`または`manual`: `n_flagella`は0以上の整数

## QC

Phase 3はclipごとに少なくとも次を記録する．

* `qc.status`: `pass`，`review`，`fail`
* `qc.exclusion_reason`: 除外理由．通過時は`null`

実動画では`qc.detection_confidence_min`，`qc.tracking_gap_count`，`qc.notes`を任意で保持できる．

擬似動画のGT bbox，render availability，Phase 2 shape QCも，共通の`qc.status`と`qc.exclusion_reason`へ反映する．

## Sources Of Truth

* Machine-readable schema: `schemas/phase3_clip_metadata.schema.json`
* Minimal fixture: `examples/phase3/clip_metadata_minimal.json`
* Contract test: `tests/test_phase3_clip_metadata_schema.py`
* Current state: `docs/phase3/phase3_current.md`
* Adopted decisions: `docs/phase3/phase3_tasks.md`
