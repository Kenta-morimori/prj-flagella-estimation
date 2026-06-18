# Phase 2.8: べん毛数分析用特徴量定義

この文書は，Issue #65「べん毛数分析用の特徴量の定義」の完了成果物である。

親Issue #71「べん毛数の違いによる遊泳挙動の分析」では，RUN固定条件の3Dシミュレーション結果を使い，べん毛数差が遊泳挙動へ与える影響をデータセットとして評価する。本Issueでは，後続の複数run実行，dataset化，分布可視化に先立ち，特徴量カテゴリと利用方針を固定する。

## 対象と前提

初期分析は以下を基本条件とする。

* `n_flagella = 1, 2, 3, 6`
* RUN固定条件
* `motor.enable_switching=false`
* `motor.force_distribution=material_twist_local_couple`
* `time.dt_star=1.0e-4`
* `flagella.initial_helix_axis_from_rear_deg=0`
* 初期評価時間は `duration_s=0.5` を基本にする

1 sample は，1つの `n_flagella`，1つの `seed`，1つの simulation run の組として定義する。

## Feature Registry

以下の YAML は，特徴量カテゴリと代表変数名を示す簡易 registry である。厳密な算出式，単位，追加派生特徴量は，後続のdataset構築タスクで定義する。

```yaml
feature_categories:
  metadata:
    target: sample_and_simulation_settings
    ml_candidate: false
    variables:
      - sample_id
      - dataset_id
      - n_flagella
      - seed
      - duration_s
      - dt_star
      - torque_Nm
      - force_distribution
      - condition_tag

  quality:
    target: qc
    ml_candidate: false
    variables:
      - quality_class
      - shape_pass
      - relaxed_pass
      - use_for_analysis
      - use_for_ml_candidate
      - review_required
      - valid_duration_s

  cell_translation:
    target: cell
    ml_candidate: true
    variables:
      - cell_displacement
      - cell_path_length
      - cell_mean_speed
      - cell_speed_std
      - cell_speed_cv
      - cell_straightness

  cell_orientation:
    target: cell
    ml_candidate: true
    variables:
      - cell_axis_angle_change
      - cell_axis_angle_std
      - cell_angular_velocity_mean
      - cell_angular_velocity_std
      - cell_angular_velocity_rms
      - cell_wobble

  flagella_axis:
    target: flagella
    ml_candidate: true
    variables:
      - flagella_axis_alignment
      - flagella_axis_spread
      - flagella_axis_pair_angle_mean
      - flagella_axis_pair_angle_max
      - flagella_axis_rear_alignment

  cell_flagella_relation:
    target: relation
    ml_candidate: true
    variables:
      - cell_flagella_axis_angle
      - cell_flagella_axis_angle_std
      - cell_flagella_axis_stability

  diagnostics:
    target: qc
    ml_candidate: false
    variables:
      - hook_drift
      - hook_wrapped
      - flyaway
      - abnormal_rotation
      - first_fail_category
```

## 利用方針

特徴量は，まず分析用特徴量として使い，その一部だけをML入力候補にする。

`quality` と `diagnostics` は，QC，外れ値確認，代表サンプル選定，異常条件の記録に使う。これらは原則としてML入力候補には含めない。ただし，後続タスクでデータ除外や層別集計に使うことは許容する。

`n_flagella=1` では，べん毛軸同士のpair angleやspreadなど，複数べん毛を前提にする特徴量は定義できない。このような特徴量は `NaN` として出力する。summaryやplotでは，`NaN` の除外数と定義不能である理由を記録する。

`use_for_analysis` と `use_for_ml_candidate` は分けて扱う。hook wrapped や relaxed pass を含むsampleは分析対象に残せる場合があるが，ML候補として使うかは後続の品質基準で判断する。

## 出力仕様

出力先は，タイムスタンプだけに依存せず，明示的な `dataset_id` を持つ以下の形式を正とする。

```text
outputs/analysis/flagella_count_behavior/datasets/<dataset_id>/summary.csv
outputs/analysis/flagella_count_behavior/datasets/<dataset_id>/timeseries/<sample_id>.csv
```

`summary.csv` は 1 sample = 1 row とし，`metadata`，`quality`，集計特徴量，診断特徴量を含める。`timeseries/<sample_id>.csv` は，sample単位の時系列を保持し，後続の特徴量再計算や外れ値確認に使う。

後続の #71 実装では，複数条件実行，dataset構築，分布可視化をこの registry と出力仕様に従って追加する。
