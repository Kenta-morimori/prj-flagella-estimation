# Phase 2 to Phase 3 handoff

作成日: 2026-07-21

## Summary

Issue #125 では，Phase 2 の3D RUN固定 dataset v1 を Phase 3 へ渡す baseline と，Phase 3 / 4 で詰める後続タスクの境界を整理する。

この文書は移行計画の正本であり，2D識別性評価，個体clip schema実装，長時間run，追加dataset生成は完了範囲に含めない。

## Handoff baseline

Phase 3 MVP へ渡す baseline は次とする。

| item | value |
| --- | --- |
| source dataset | Phase 2 analysis dataset v1 |
| config | `conf/phase2_multi_run/flagella_count_behavior_v1.yaml` |
| raw run root | `outputs/phase2_multi_run/flagella_count_behavior_v1` |
| analysis dataset | `outputs/phase2_analysis/flagella_count_behavior/datasets/v1` |
| model_id | `flag_spring2p25_body2p5_candidate` |
| state | RUN fixed, `motor.enable_switching=false` |
| candidate labels | `n_flagella=1,2,3` |
| candidate samples | 27 strict-pass samples |
| split group key | `run_id` or a stable ID derived from one simulation condition |

`n_flagella=1,2,3` は全27 sampleが全時系列 strict pass し，2026-07-16 の user visual review で body deformation，helical shape preservation，swimming-like propulsion に training candidate から外す問題はないと判断された。

## Baseline exclusions

次は Phase 3 MVP の blocker にしない。

| condition | handoff status | follow-up |
| --- | --- | --- |
| `n_flagella=4` | diagnostic-only | Issue #124 |
| `n_flagella>=5` | excluded | Issue #124 |
| TUMBLE / motor switching | baseline外 | Issue #69 |
| Brownian項 | baseline外 | Issue #15 |
| torque variation | augmentationとして未採用 | Issue #128 |
| model equation / stiffness profile変更 | dataset version変更候補 | Issue #128 |
| 実顕微鏡入力条件 | Phase 3側で整理 | Issue #8 / #9 |
| 2D描画条件 | Phase 2→3側で整理 | Issue #17 |

`n_flagella=4` は 9 sample中6 sample が strict pass したが，3 sample に `flag` failure があり，replayでも非等速回転が残った。そのため training candidate へは入れず，物理改善は #124 へ分離する。

## Data units

Phase 3 / 4 では，次の単位を区別する。

| term | definition | notes |
| --- | --- | --- |
| `run` | 1つの3D simulation conditionの実行結果 | `n_flagella`, `attach_seed`, `phase_seed`, model conditionで一意化する |
| `source_video` | crop前の入力動画 | 擬似動画では Phase 2 render，実動画では顕微鏡動画 |
| `track` | 1個体の連続観測 | 擬似単一菌体では1 runに1 trackを基本とする |
| `clip` | Phase 4 の学習入力候補となる固定長 frame sequence | 同一 track から複数生成しても独立run数には数えない |
| `frame` | 1時刻の画像と対応metadata | bbox，時刻，projection / normalization情報を持つ |

独立データ数の基本単位は `run` または実動画の独立 `track` とする。同一 `run` / `track` 由来の複数 `clip` は train / validation / test をまたがせない。

## ID responsibilities

| id | responsibility |
| --- | --- |
| `model_id` | 物理モデル・主要stiffness・force distributionなど，dataset解釈に影響する条件を表す |
| `dataset_id` | analysis dataset の論理版を表す。v1 は今回の handoff baseline |
| `run_id` | 1つのsimulation condition / 実track sourceを識別する |
| `render_id` | 同じ run から作られた描画条件を識別する。noise / blur / projection / fps などを含める |
| `source_video_id` | crop前動画を識別する。実動画では取得動画，擬似動画では render 出力 |
| `track_id` | source video内の個体trackを識別する |
| `clip_id` | trackから切り出した学習入力候補を識別する |
| `frame_id` | source videoまたはclip内の1 frameを識別する |

Split leakage を防ぐため，Phase 4 の group key は少なくとも `run_id` または `track_id` を含める。`render_id` や `clip_id` だけで split してはいけない。

## Input paths

Phase 3 は入力経路を分けるが，共通の clip / metadata schema へ収束させる。

| input | path | detection |
| --- | --- | --- |
| 実顕微鏡動画 | source video -> detection -> tracking -> crop -> normalization -> clip | required |
| Phase 2 擬似動画 | source video + GT bbox / track -> crop -> normalization -> clip | skippable |
| 将来の複数菌体擬似動画 | GT passthrough または detection -> tracking -> crop -> normalization -> clip | selectable |

擬似動画で detection を省略しても，Phase 3 の出力 schema は実動画経路と同じにする。実動画に存在しない `n_flagella` label は unknown とし，擬似動画の GT label と混同しない。

## Immediate follow-up issues

| issue | role | dependency |
| --- | --- | --- |
| #126 | 2D動画上でも `n_flagella=1,2,3` の差が残るか確認する | #125 / #119 |
| #127 | 実動画・擬似動画の共通clip schemaを決める | #125 |
| #128 | augmentation / domain variation / dataset version変更の規則を決める | #125 / #119 |
| #129 | clip時間長と必要独立run数を決める | #126 / #127 |
| #6 | Phase 3 clip生成パイプラインを実装する | #125 / #127 |

長期実行は #126 と #129 で扱う。#125 では長期simulation，full 2D render，大量clip生成は実行しない。

## User decisions carried forward

現時点では次の判断で進める。

1. Phase 3 MVP は v1 `n_flagella=1,2,3` を baseline とする。
2. `n_flagella>=4` の改善は Phase 3 MVP を block しない。
3. まず擬似動画の GT passthrough 経路を優先し，実動画 detection / tracking は同じ出力schemaへ後続実装する。
4. 長期実行はユーザーが実施し，Codex はコマンド，出力先，評価観点，smoke checkを用意する。

これらを変更する場合は，#125ではなく #126 / #127 / #128 / #129 または #6 の開始時に更新する。
