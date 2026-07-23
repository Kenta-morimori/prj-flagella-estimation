# Phase 3.2 clip duration and independent run count design

この文書は Issue #129 の設計メモである。Phase 3 / 4 MVP の標準 clip 時間長を `0.5 s` とし，評価単位，split 規則，軽量 fixture，重い learning curve 実行の渡し方を固定する。

## Scope

ここで決める:

- `clip duration` / `source duration` / `simulation run` / `track` / `window` の用語
- #127 schema の `track.group_key` を使う grouped split 規則
- `0.25 s` / `0.5 s` / `1.0 s`，non-overlap / overlap の評価設計
- 独立データ数と反復観測数の数え方
- 重い learning curve 実行をユーザーへ渡す条件

ここで決めない:

- `n_flagella` ごとの必要 run 数
- 実動画 detection / tracking の採用手法
- Phase 2 simulation の物理条件変更

## Fixed Decision

Phase 3 / 4 MVP の標準 clip duration は `0.5 s` とする。

- `0.5 s` は Phase 2 v1 の `1.0 s` source run から non-overlap で2 windowを切れるため，現行出力を使った pilot に適している。
- `0.25 s` は短い実trackへの耐性確認，`1.0 s` は現行 full-run baseline との比較として残す。
- 採用後も `0.25 s` / `1.0 s` は ablation / robustness 条件であり，MVP dataset の default ではない。
- overlap window は inference-time 安定性や反復観測の評価には使うが，独立 sample 数には数えない。

## Terms

| term | definition | independent sample? |
| --- | --- | --- |
| `simulation run` | Phase 2 の1条件・1 seed 組み合わせから得る raw simulation。例: `as000__ps000__nf01` | yes for pseudo data |
| `source video` | render 済み動画または実顕微鏡動画。1つの `simulation run` から複数 `render_id` の source video が作られうる | no by itself |
| `source duration` | source video 全体の時間長。複数 window を切るため，評価対象の clip duration 以上にする | no |
| `track` | 1個体の時系列。Phase 2 GT passthrough では通常 `gt_body` の1 track，実動画では tracker 出力の1 track | yes for real data; pseudo では run と合わせて数える |
| `window` | track から切り出す frame interval。non-overlap または overlap のどちらか | no |
| `clip duration` | Phase 4 入力候補 clip の時間長。MVP default は `0.5 s`，比較候補は `0.25 s`, `1.0 s` | no |
| `clip` | window から crop / normalization 後に出力される frame sequence | no |

`clip` 数や frame 数を独立データ数として数えない。独立データ数は pseudo では `simulation run`，実動画では `source_video_id` と `track_id` の組み合わせを基本単位にする。

## Grouped Split Rule

Phase 4 split は #127 schema の `track.group_key` を唯一の group 正本として使う。train / validation / test のどれかに割り当てた `group_key` は，他 split へ出してはいけない。

Pseudo data の推奨 `group_key`:

```text
phase2:<dataset_version>:<run_id>
```

実動画の推奨 `group_key`:

```text
real:<source_video_id>:<track_id>
```

同一 `run_id` から作られた `render_id` 違い，同一 track から作った overlap window，正規化設定だけが異なる clip は同じ `group_key` に置く。`clip_id`, `render_id`, `window_policy` だけで split してはいけない。

## Window Policies

| policy | stride | purpose | independent count |
| --- | --- | --- | --- |
| `non_overlap` | `clip_duration_s` | 評価の基本。track 内相関を抑えた window 切り出し | 追加しない |
| `overlap` | `clip_duration_s / 2` for the initial comparison | inference-time の頑健性や反復観測として比較 | 追加しない |
| `full_run` | source 全体 | 後方互換 smoke / 1.0 s baseline | 追加しない |

overlap window は，情報量や安定性の比較には使えるが，learning curve の x 軸には足さない。metadata の `clip.window_policy` は #127 schema に合わせて `overlap` とし，50% stride は config / manifest 側の window stride 設定で表す。learning curve の x 軸は `n_flagella` ごとの unique `group_key` 数にする。

## Initial Evaluation Matrix

| clip duration | non-overlap | overlap | required source duration |
| ---: | --- | --- | --- |
| `0.25 s` | ablation | `50%` overlap comparison | `>= 1.0 s` pilot; longer run if variance is high |
| `0.5 s` | MVP default | `50%` overlap comparison | `>= 1.0 s` pilot; longer run if non-overlap count is insufficient |
| `1.0 s` | current full-run baseline | no overlap in current `1.0 s` source | `>= 1.0 s`; longer source needed for multi-window |

初期 metric:

- class-balanced accuracy / macro F1
- group-level bootstrap confidence interval
- `n_flagella` ごとの confusion matrix
- clip 内 frame 由来特徴量の安定性
- 同一 group 内 overlap / non-overlap 予測の分散
- `qc.status=pass` 比率と exclusion reason 分布

## Lightweight Implementation Boundary

#6 の最小 pipeline に渡すため，まず次の小さな fixture / library test で固定する。

- #127 の `examples/phase3/clip_metadata_minimal.json` から `clip.duration_s` と `clip.window_policy` のバリエーションを生成できる。
- Phase 3 pipeline default は `clip.duration_s=0.5`, `clip.window_policy=non_overlap` とする。
- window 切り出し関数は frame index と fps から `source_frame_start`, `source_frame_end`, `t_start_s`, `t_end_s` を決める。
- split helper は `track.group_key` 単位で split し，同じ group が複数 split に入らないことを test する。
- 実動画の detection / tracking は fixture では mock track を使い，実AVI解析や detector 採用判断を含めない。

## User-Run Required: Grouped Learning Curve

次の段階では，Phase 2 dataset v1 の既存 pseudo dataset から clip duration 比較用 dataset を作り，grouped learning curve を評価する必要がある。これは重い実行・採用判断を含むため，Codex は勝手に実行しない。

Exact command draft:

```bash
uv run python scripts/03_phase3/build_clip_dataset.py \
  config=conf/phase3/clip_duration_pilot.yaml \
  input_dataset=outputs/phase2_analysis/flagella_count_behavior/datasets/v1 \
  output_dir=outputs/YYYY-MM-DD/HHMMSS/phase3_clip_duration_pilot \
  clip_durations_s='[0.25,0.5,1.0]' \
  window_policies='[non_overlap,overlap]' \
  window.overlap_stride_fraction=0.5 \
  split.group_key_field=track.group_key

uv run python scripts/04_train_flagella_count/evaluate_learning_curve.py \
  dataset_dir=outputs/YYYY-MM-DD/HHMMSS/phase3_clip_duration_pilot \
  output_dir=outputs/YYYY-MM-DD/HHMMSS/phase4_grouped_learning_curve_clip_duration \
  group_key_field=track.group_key \
  class_field=labels.n_flagella
```

Expected output files:

```text
outputs/YYYY-MM-DD/HHMMSS/phase3_clip_duration_pilot/manifest.json
outputs/YYYY-MM-DD/HHMMSS/phase3_clip_duration_pilot/clip_metadata.jsonl
outputs/YYYY-MM-DD/HHMMSS/phase3_clip_duration_pilot/split_summary.csv
outputs/YYYY-MM-DD/HHMMSS/phase4_grouped_learning_curve_clip_duration/learning_curve.csv
outputs/YYYY-MM-DD/HHMMSS/phase4_grouped_learning_curve_clip_duration/metrics_by_duration.csv
outputs/YYYY-MM-DD/HHMMSS/phase4_grouped_learning_curve_clip_duration/confusion_by_duration.csv
outputs/YYYY-MM-DD/HHMMSS/phase4_grouped_learning_curve_clip_duration/manifest.json
```

Decision points:

- `0.5 s` default が `0.25 s` と比べて class confusion を抑え，`1.0 s` と比べて source duration 要求を抑えられるか。
- `0.25 s` が短い実trackへの耐性を上げる一方で class confusion を増やしすぎないか。
- `1.0 s` が性能を上げる一方で独立 group 数不足や source duration 要求を大きくしすぎないか。
- overlap window が評価安定化に効くか。ただし独立 sample 数として数えない前提を守れるか。
- `n_flagella=1,2,3` の各 class で追加 run がどの程度必要か。
- 実動画 track 長の想定と矛盾しない clip duration か。

上記 command は現時点では設計上の draft であり，#6 / Phase 4 側の CLI 実装後に実コマンドとして確定する。
