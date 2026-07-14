# Phase 2.8 diagnostic dataset v0 report

作成日: 2026-07-14

対象:

- Issue #117: 現状モデルの diagnostic dataset v0 統計レポート
- Issue #118: `model_id` ベースの analysis dataset config / `dataset_id` 命名規則

## Positioning

diagnostic dataset v0 は Phase3/4 学習用 dataset の凍結ではない。現状モデルで RUN 固定べん毛本数差が特徴量へ出るか，また `n>=4` の破綻境界がどこにあるかを確認する diagnostic baseline である。

統計 source は既存 historical alias dataset:

- config: `conf/phase2_multi_run/flagella_count_behavior_diagnostic.yaml`
- dataset_id: `fc_nf1_2_3_6_as3_ps3_torque2p0_dur1p0`
- dataset output: `outputs/phase2_analysis/flagella_count_behavior/datasets/fc_nf1_2_3_6_as3_ps3_torque2p0_dur1p0`
- raw run root: `outputs/phase2_multi_run/flagella_count_behavior_diagnostic`

Issue #118 以降の canonical v0 entrypoint は次にする:

- config: `conf/phase2_multi_run/flagella_count_behavior_runfixed_rtseg_fp1p25_torque2p0_v0.yaml`
- dataset_id: `fc_runfixed_rtseg_fp1p25_torque2p0_v0_nf1_2_3_6_as3_ps3_dur1p0`
- dataset output: `outputs/phase2_analysis/flagella_count_behavior/datasets/fc_runfixed_rtseg_fp1p25_torque2p0_v0_nf1_2_3_6_as3_ps3_dur1p0`
- raw run root: `outputs/phase2_multi_run/flagella_count_behavior_runfixed_rtseg_fp1p25_torque2p0_v0`

今回の統計レポートでは canonical v0 output は再生成していない。既存 output と docs/tests の互換性を保つため，historical alias は残す。

## Conditions

v0 の条件:

- `n_flagella=[1,2,3,6]`
- `attach_seed=[0,1,2]`
- `phase_seed=[0,1,2]`
- `duration_s=1.0`
- `time.dt_star=1.0e-4`
- `motor.torque_Nm=2.0e-20`
- `motor.enable_switching=false`
- `motor.force_distribution=root_torque_segment_couples`
- `flagella.placement_mode=seeded_surface`
- `flagella.initial_phase_mode=seeded`
- `flagella.initial_helix_axis_from_rear_deg=0`

`model_id` は `runfixed_rtseg_fp1p25_torque2p0` とする。`fp1p25` は #103 後の basal freedom default である `motor.local_attach_frame_position_scale=1.25` を指す。`n_flagella`，seed 数，duration は model ではなく dataset 条件として `dataset_id` の後半に置く。

## QC Summary

| `n_flagella` | sample count | `quality_class` | `first_fail_category` | `use_for_ml_candidate` | v0 interpretation |
| --- | ---: | --- | --- | --- | --- |
| 1 | 9 | all `strict_pass` | all `none` | 9/9 true | diagnostic baseline / ML candidate range |
| 2 | 9 | all `strict_pass` | all `none` | 9/9 true | diagnostic baseline / ML candidate range |
| 3 | 9 | all `strict_pass` | all `none` | 9/9 true | diagnostic baseline / ML candidate range |
| 6 | 9 | all `fail` | all `flag` | 0/9 true | diagnostic-only failure condition |

`n=1,2,3` は v0 の範囲では transient fail を示していない。今後 v1 や追加条件で `relaxed_pass` だが `strict_pass` でない transient sample が出た場合は，Phase3/4 training candidate へ直行させず，analysis-only または user visual review required として扱う。

## Feature Separation

代表特徴量の mean +/- sample std:

| feature | n=1 | n=2 | n=3 | n=6 |
| --- | ---: | ---: | ---: | ---: |
| `cell_mean_speed` | 0.1896 +/- 0.0134 | 0.3767 +/- 0.1186 | 0.4151 +/- 0.0772 | 0.6464 +/- 0.0802 |
| `cell_straightness` | 0.3225 +/- 0.0281 | 0.1775 +/- 0.0085 | 0.1559 +/- 0.0490 | 0.1205 +/- 0.0222 |
| `cell_angular_velocity_rms` | 1.3013 +/- 0.1751 | 2.3178 +/- 0.1166 | 5.7270 +/- 2.2220 | 10.6868 +/- 1.5020 |
| `hook_drift` | 0.0216 +/- 0.0132 | 0.1786 +/- 0.0126 | 0.1206 +/- 0.0503 | 0.1316 +/- 0.0241 |
| `flagella_axis_alignment` | NA | 0.9435 +/- 0.0008 | 0.9570 +/- 0.0179 | 0.9653 +/- 0.0163 |

`n=1,2,3` の特徴分離は，とくに `cell_mean_speed`，`cell_straightness`，`cell_angular_velocity_rms`，`hook_drift` で見える。`flagella_axis_alignment` は単一べん毛では pairwise alignment が定義されないため `n=1` は NA とする。

feature screening では，上位の連続特徴として次が分離候補だった:

| feature | range / within std |
| --- | ---: |
| `hook_drift` | 6.565 |
| `cell_straightness` | 6.196 |
| `cell_axis_angle_std` | 5.725 |
| `cell_angular_velocity_rms` | 5.602 |
| `cell_axis_angle_change` | 5.349 |
| `cell_angular_velocity_mean` | 5.349 |
| `cell_path_length` | 5.041 |
| `cell_mean_speed` | 3.431 |

ただし `n=6` は fail sample なので，これらの分離は `n=1,2,3` の diagnostic baseline として読む。`n=6` の高速化や角速度増加は，training label として信頼するのではなく破綻診断値として扱う。

## n>=4 Interpretation

Issue #113 / PR #114 の seed 固定診断では，`attach_seed=0`, `phase_seed=0`, `n_flagella=4,5,6` がすべて `flag` first fail となった。body diagnostics 対応後の再実行では，全条件で `body_shape_pass=false`, `body_fail_category=body_spring` となり，ユーザー定性評価でも body 伸長が見えると判断された。

| `n_flagella` | `first_fail_t_s` | `max_flag_bond_rel_err` | `body_fail_category` | `body_spring_max_stretch_ratio` |
| --- | ---: | ---: | --- | ---: |
| 4 | 0.3273 | 1.8486 | `body_spring` | 1.4972 |
| 5 | 0.3168 | 1.5037 | `body_spring` | 1.3724 |
| 6 | 0.3505 | 1.4889 | `body_spring` | 1.3863 |

このため，現状モデルでは `n>=4` を Phase3/4 training candidate から外し，diagnostic-only とする。`n>=4` の復帰判断は #115 のモデル修正と #119 の改善モデル dataset v1 再生成後に行う。

## Naming Rule

analysis dataset config は次の規則で命名する:

```text
conf/phase2_multi_run/flagella_count_behavior_<model_id>_<dataset_version>.yaml
```

analysis dataset id は次の規則で命名する:

```text
fc_<model_id>_<dataset_version>_nf<flagella-list>_as<attach-seed-count>_ps<phase-seed-count>_dur<duration>
```

`model_id` に含める主要条件:

- motor state: `runfixed`
- torque force distribution: `rtseg`
- basal freedom representative: `fp1p25`
- torque: `torque2p0`

`model_id` に含めない dataset 条件:

- `n_flagella`
- `attach_seed`
- `phase_seed`
- `duration_s`

v1 では改善モデルの主要変更を `model_id` に反映し，`dataset_version=v1` を使う。v1 の `n_flagella` 範囲は #115/#119 の結果で決める。

## Next Actions

- #115: `n>=4` の flag bond 過伸長と body_spring failure をモデル側で改善する。
- #116: #115 の候補が出た後，少数条件 sweep / heatmap 方針を決める。
- #119: 改善モデルで analysis dataset v1 を再生成し，v0 の `n=1,2,3` baseline と比較する。
