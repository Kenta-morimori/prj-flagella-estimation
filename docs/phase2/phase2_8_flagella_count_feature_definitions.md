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

Feature Registry の正本は `conf/phase2_analysis/flagella_count_behavior_features.yaml` とする。

この YAML は，特徴量カテゴリ，対象領域，ML入力候補フラグ，代表変数名，短い説明コメントを定義する。厳密な算出式，単位，追加派生特徴量，`summary.csv` / `timeseries` の完全schemaは，後続のdataset構築タスクで定義する。

正本YAMLでは，各カテゴリを以下の形式で記録する。

```yaml
feature_categories:
  category_name:
    target: cell | flagella | relation | qc | sample_and_simulation_settings
    ml_candidate: true | false
    variables:
      - representative_variable_name  # 短い説明
```

現在のカテゴリは，`metadata`, `quality`, `cell_translation`, `cell_orientation`, `flagella_axis`, `cell_flagella_relation`, `diagnostics` である。各カテゴリの代表変数名は，`conf/phase2_analysis/flagella_count_behavior_features.yaml` を参照する。

## 利用方針

特徴量は，まず分析用特徴量として使い，その一部だけをML入力候補にする。

`quality` と `diagnostics` は，QC，外れ値確認，代表サンプル選定，異常条件の記録に使う。これらは原則としてML入力候補には含めない。ただし，後続タスクでデータ除外や層別集計に使うことは許容する。

`n_flagella=1` では，べん毛軸同士のpair angleやspreadなど，複数べん毛を前提にする特徴量は定義できない。このような特徴量は `NaN` として出力する。summaryやplotでは，`NaN` の除外数と定義不能である理由を記録する。

`use_for_analysis` と `use_for_ml_candidate` は分けて扱う。hook wrapped や relaxed pass を含むsampleは分析対象に残せる場合があるが，ML候補として使うかは後続の品質基準で判断する。

## 出力仕様

本Issueでは，dataset保存先の具体パスは決めない。複数runを束ねるdatasetの保存先，raw simulation runとの対応，`run.log` / `manifest.json` の配置は，親Issue #71 のdataset構築タスクで設計する。

後続タスクでは，少なくとも以下を満たす必要がある。

* `summary.csv` は 1 sample = 1 row とし，`metadata`，`quality`，集計特徴量，診断特徴量を含める。
* `timeseries/<sample_id>.csv` は，sample単位の時系列を保持し，後続の特徴量再計算や外れ値確認に使う。
* `dataset_id` は，timestampとは別の論理IDとして持つ。
* `sample_id` は，1つの `n_flagella`，1つの `seed`，1つの simulation run の組を識別できる形にする。
* dataset生成条件，raw run paths，config paths，overrides，seed，Git commit情報を再現性metadataとして記録する。

後続の #71 実装では，複数条件実行，dataset構築，分布可視化を `conf/phase2_analysis/flagella_count_behavior_features.yaml` と出力仕様に従って追加する。
