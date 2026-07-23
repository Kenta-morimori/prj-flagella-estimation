# Phase 2.8 2D projection separability

作成日: 2026-07-22

## Summary

Issue #126 では，Phase 2 dataset v1 の RUN 固定 `n_flagella=1,2,3` について，3Dで見えていた本数差が2D投影後の運動特徴量にも残るかを初期確認した。

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

初期実装では，2D動画そのものの画像特徴ではなく，`trajectory.csv` から2D投影で観測できる運動特徴量を抽出する。dataset v1 の2D rendererは `render.center_body_in_2d=true` で各frameを body center に合わせるため，body center の raw XY translation は `projection.mp4` では見えない。したがって baseline では renderer camera frame に合わせ，body center translation 由来の特徴量はゼロとして扱い，body axis のXY投影角度特徴を主な情報として使う。

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

`n_flagella=1,2,3` の27 sampleでは，renderer camera frame に合わせた grouped nearest-centroid baseline の accuracy は `19/27 = 0.704` だった。

誤分類:

| sample_id | actual | predicted |
| --- | ---: | ---: |
| `nf01_as001_ps000` | 1 | 2 |
| `nf02_as001_ps000` | 2 | 1 |
| `nf03_as001_ps000` | 3 | 2 |
| `nf01_as001_ps001` | 1 | 2 |
| `nf02_as001_ps002` | 2 | 1 |
| `nf03_as001_ps002` | 3 | 2 |
| `nf02_as002_ps000` | 2 | 1 |
| `nf02_as002_ps002` | 2 | 1 |

主な2D特徴量の平均:

| feature | n=1 | n=2 | n=3 |
| --- | ---: | ---: | ---: |
| `xy_displacement_um` | 0.0000 | 0.0000 | 0.0000 |
| `xy_path_length_um` | 0.0000 | 0.0000 | 0.0000 |
| `xy_mean_speed_um_s` | 0.0000 | 0.0000 | 0.0000 |
| `xy_body_axis_angle_std_deg` | 9.7087 | 13.2245 | 13.0673 |
| `xy_body_axis_angular_velocity_rms_rad_s` | 0.9550 | 1.7708 | 5.2037 |

この初期結果では，body center を固定した2D camera frameでも body axis の向き変化には `n_flagella` による差が一部残る。ただし raw translation を除くと識別性は `26/27` から `19/27` へ下がるため，実際の2D pseudo-microscopy frameから抽出できる画像特徴や短いclip窓での識別性を #127 / #129 で評価する必要がある。

## Observation Scale

今回の解析単位は dataset v1 の各 raw condition 全体であり，source run duration は `1.0 s` である。同一 run 由来の frame / clip は独立sampleとして数えず，`attach_seed` / `phase_seed` の組を group として扱う。

Issue #126 の完了判断は「1.0 s run 全体の body-centered 2D camera frame で，本数差が初期特徴量にどの程度残るか」の確認までとする。Phase 4 の入力clip長は別設計値であり，`0.25 s` / `0.5 s` / `1.0 s` の短時間窓比較と必要な独立run数は #129 へ引き渡す。

## Handoff

次の確認対象:

- `projection.mp4` / frame sequence から同じ傾向を取り出せるか。
- 0.25 s / 0.5 s / 1.0 s clip で2D特徴量がどの程度安定するか。
- grouped split を維持した baseline を #129 の clip時間長評価へ接続する。
- Phase 3共通 clip / metadata schema では，`run_id` または `source_video_id` / `track_id` を group key として保持する。
