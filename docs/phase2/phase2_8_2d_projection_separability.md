# Phase 2.8 2D projection separability

作成日: 2026-07-22

## Summary

Issue #126 では，Phase 2 dataset v1 の RUN 固定 `n_flagella=1,2,3` について，3Dで見えていた本数差がXY投影後の運動特徴量にも残るかを初期確認した。

- source dataset: `outputs/phase2_analysis/flagella_count_behavior/datasets/v1`
- target samples: `use_for_ml_candidate=True` かつ `n_flagella=1,2,3`
- sample count: 27
- input: 各 raw condition の `trajectory.csv`
- output: `outputs/phase2_analysis/flagella_count_behavior/datasets/v1/analysis/projection_2d`

## Command

```bash
uv run python scripts/02_phase2_analysis/analyze_2d_separability.py \
  --dataset-dir outputs/phase2_analysis/flagella_count_behavior/datasets/v1 \
  --output-dir outputs/phase2_analysis/flagella_count_behavior/datasets/v1/analysis/projection_2d \
  --overwrite
```

出力:

- `features_2d.csv`: sampleごとの2D特徴量
- `feature_summary_by_n_flagella.csv`: `n_flagella` 別の要約統計
- `grouped_nearest_centroid_baseline.csv`: `attach_seed` / `phase_seed` group leave-out の最近傍centroid baseline
- `manifest.json`: 入力・出力・対象sample数・baseline accuracy

## Feature Scope

初期実装では，2D動画そのものの画像特徴ではなく，`trajectory.csv` からXY投影で観測できる運動特徴量を抽出する。

- `xy_displacement_um`
- `xy_path_length_um`
- `xy_mean_speed_um_s`
- `xy_speed_std_um_s`
- `xy_speed_cv`
- `xy_straightness`
- `xy_range_x_um`
- `xy_range_y_um`
- `xy_extent_um`
- `xy_body_axis_angle_change_deg`
- `xy_body_axis_angle_std_deg`
- `xy_body_axis_angular_velocity_rms_rad_s`
- `xy_velocity_heading_change_deg`
- `xy_velocity_heading_std_deg`

同一 simulation run 由来の frame / clip を独立sampleとして数えない方針に合わせ，簡易 baseline は `attach_seed` と `phase_seed` の組を group として leave-one-group-out 評価する。

## Initial Result

`n_flagella=1,2,3` の27 sampleでは，grouped nearest-centroid baseline の accuracy は `26/27 = 0.963` だった。

唯一の誤分類:

| sample_id | actual | predicted |
| --- | ---: | ---: |
| `nf02_as001_ps002` | 2 | 3 |

主な2D特徴量の平均:

| feature | n=1 | n=2 | n=3 |
| --- | ---: | ---: | ---: |
| `xy_displacement_um` | 0.1507 | 0.2550 | 0.3381 |
| `xy_path_length_um` | 0.4660 | 1.7896 | 2.8266 |
| `xy_mean_speed_um_s` | 0.1507 | 0.2550 | 0.3381 |
| `xy_straightness` | 0.3259 | 0.1390 | 0.1288 |

この初期結果では，XY投影後も `n_flagella` による運動分布差は残っている。ただし，まだ `trajectory.csv` ベースの幾何特徴であり，実際の2D pseudo-microscopy frameから抽出できる画像特徴や短いclip窓での識別性は未評価である。

## Observation Scale

今回の解析単位は dataset v1 の各 raw condition 全体であり，source run duration は `1.0 s` である。同一 run 由来の frame / clip は独立sampleとして数えず，`attach_seed` / `phase_seed` の組を group として扱う。

Issue #126 の完了判断は「1.0 s run 全体のXY投影特徴量で本数差が残るか」の初期確認までとする。Phase 4 の入力clip長は別設計値であり，`0.25 s` / `0.5 s` / `1.0 s` の短時間窓比較と必要な独立run数は #129 へ引き渡す。

## Handoff

次の確認対象:

- `projection.mp4` / frame sequence から同じ傾向を取り出せるか。
- 0.25 s / 0.5 s / 1.0 s clip で2D特徴量がどの程度安定するか。
- grouped split を維持した baseline を #129 の clip時間長評価へ接続する。
- Phase 3共通 clip / metadata schema では，`run_id` または `source_video_id` / `track_id` を group key として保持する。
