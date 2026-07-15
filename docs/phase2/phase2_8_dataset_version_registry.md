# Phase 2.8 dataset version registry

作成日: 2026-07-14

この文書は，Phase 2.8 の flagella-count behavior analysis dataset について，dataset version と入出力対応を一覧で管理する registry である。

config と manifest が実行条件の source of truth であり，この registry は人間が version の意図，canonical config，historical output，後続 task との関係をすばやく確認するための索引である。

## Version policy

### Fields

| field | meaning |
| --- | --- |
| `dataset_version` | 比較単位としての dataset 論理版。主系列は `v0`, `v1`, `v2` とする。 |
| `dataset_revision` | 同じ `dataset_version` の再生成，小修正，manifest/docs 整理を区別する補助版。`r0`, `r1` とする。 |
| `dataset_scope` | dataset 条件の短い要約。path ではなく registry / config / manifest で読む。 |
| `model_revision` | 物理モデル・主要パラメータ条件の短い論理ID。詳細値は config / manifest を読む。 |
| `version_role` | `diagnostic_baseline`, `diagnostic_sweep`, `training_candidate`, `planned` などの用途。 |
| `training_candidate` | Phase3/4 学習候補として扱うか。diagnostic baseline は `false`。 |

### Version bump rules

`dataset_version` を上げる変更:

- scale / basal freedom / force distribution / torque transmission など，モデル解釈が変わる変更
- `n>=4` failure 改善モデルの採用
- Phase3/4 training candidate の本数範囲が変わる変更
- canonical dataset として v0 と直接比較する改善モデル dataset

`dataset_version` を上げない変更:

- torque 数値だけの diagnostic sweep
- `duration_s` 延長による安定性確認
- seed 数追加による統計補強
- 同条件再実行，manifest 修正，出力整理，docs 修正

`v0-1` / `v0.1` は使わない。必要な小修正・再生成は `dataset_revision: r1`, `r2` で管理する。path に revision を入れる必要がある場合だけ `datasets/v0_r1` のようにするが，通常は config / manifest / registry に留める。

## Registry

| dataset_version | dataset_revision | dataset_scope | status | version_role | training_candidate | model_revision | canonical config | canonical dataset_id | canonical output status | source output | interpretation |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `v0` | `r0` | `nf1_2_3_6_as3_ps3_dur1p0` | current diagnostic baseline | `diagnostic_baseline` | false | `current_v0` | `conf/phase2_multi_run/flagella_count_behavior_v0.yaml` | `v0` | `not_generated` | `outputs/phase2_analysis/flagella_count_behavior/datasets/fc_nf1_2_3_6_as3_ps3_torque2p0_dur1p0` | 現状モデルの diagnostic baseline。`n=1,2,3` は特徴分離を見る baseline，`n>=4` は diagnostic-only。 |
| `v1` | `r0` | TBD | planned | `planned` | TBD | TBD after #116 | TBD | TBD | `not_generated` | none | #115 では `n>=4` を dataset v1 へ直接戻せる条件は未確定。#116 の追加 sweep 後に #119 で再生成可否を判断する。 |

## v0 details

canonical v0 config:

- `conf/phase2_multi_run/flagella_count_behavior_v0.yaml`

historical alias:

- config: `conf/phase2_multi_run/flagella_count_behavior_diagnostic.yaml`
- dataset_id: `fc_nf1_2_3_6_as3_ps3_torque2p0_dur1p0`
- dataset output: `outputs/phase2_analysis/flagella_count_behavior/datasets/fc_nf1_2_3_6_as3_ps3_torque2p0_dur1p0`
- raw run root: `outputs/phase2_multi_run/flagella_count_behavior_diagnostic`

v0 の主要条件は config の `metadata.model_notes` と `base_overrides` を読む。summary としては次である。

- `motor.enable_switching=false`
- `motor.force_distribution=root_torque_segment_couples`
- `motor.local_attach_frame_position_scale=1.25`
- `motor.torque_Nm=2.0e-20`
- `n_flagella=[1,2,3,6]`
- `attach_seed=[0,1,2]`
- `phase_seed=[0,1,2]`
- `duration_s=1.0`

canonical output `outputs/phase2_analysis/flagella_count_behavior/datasets/v0` は，2026-07-14 時点では未生成である。既存統計レポートは historical alias output を source として読む。

## v1 handoff

v1 は #115/#116 のモデル改善結果を受けて #119 で作成する。v1 config は原則として `conf/phase2_multi_run/flagella_count_behavior_v1.yaml` とし，dataset_id は `v1` とする。#115 の `flag_spring/body` 比較では `n>=4` を直接戻せる安定条件が出ていないため，現時点では v1 を生成しない。

v1 に上げるか，v0 の `dataset_revision` として扱うかは次で判断する。

- torque / duration / seed だけの補助診断なら `dataset_version` は上げず，diagnostic sweep または `dataset_revision` で記録する。
- basal freedom scale，force distribution，failure 改善モデル，training candidate range が変わるなら `dataset_version=v1` とする。
