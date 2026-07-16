# Phase 2.8 improved-model dataset v1 report

作成日: 2026-07-16

## Summary

Issue #119 では，#116 の handoff 条件を使って analysis dataset v1 を36 sample生成した。

- config: `conf/phase2_multi_run/flagella_count_behavior_v1.yaml`
- dataset_id: `v1`
- raw run root: `outputs/phase2_multi_run/flagella_count_behavior_v1`
- dataset output: `outputs/phase2_analysis/flagella_count_behavior/datasets/v1`
- model_id: `flag_spring2p25_body2p5_candidate`
- dataset_scope: `nf1_2_3_4_as3_ps3_dur1p0`
- model condition: `stiffness_scales.flag_spring=2.25`, `stiffness_scales.body=2.5`

raw simulationは `attach_seed=0..2`, `phase_seed=0..2`, `n_flagella=1..4` の36条件すべてが `t=0.9999 s` まで完走した。dataset buildと分布分析も生成済みである。

## QC Result

dataset QCは全stepを対象にstrict / relaxed gateを評価する。途中でfailureを検出して最終stepで回復したsampleも `strict_pass` には戻さない。

| n_flagella | samples | strict_pass | fail | body pass | training interpretation |
| ---: | ---: | ---: | ---: | ---: | --- |
| 1 | 9 | 9 | 0 | 9 | training candidate |
| 2 | 9 | 9 | 0 | 9 | training candidate |
| 3 | 9 | 9 | 0 | 9 | training candidate |
| 4 | 9 | 6 | 3 | 9 | diagnostic-only |

`n=4` のfailureはすべて `flag` categoryである。

| sample_id | first_fail_t_s | max_flag_bond_rel_err | final_shape_pass_nonbody |
| --- | ---: | ---: | --- |
| `nf04_as002_ps000` | 0.3259 | 1.1506 | true |
| `nf04_as000_ps002` | 0.8971 | 1.0865 | true |
| `nf04_as002_ps002` | 0.9299 | 1.0898 | false |

最終時刻だけでは上2件をpassと誤認するため，`first_fail_t_s` を優先する。ユーザー判断により，v1のPhase3/4 training candidate範囲は `n=1,2,3` とし，`n=4` は除外する。`n=5,6` も #116 の結果に基づき対象外のままとする。

## v0 Comparison

同一の `attach_seed` / `phase_seed` を対応させた9 sampleの平均値を比較した。v0 sourceはhistorical alias dataset `fc_nf1_2_3_6_as3_ps3_torque2p0_dur1p0` である。

| n | feature | v0 mean | v1 mean | change |
| ---: | --- | ---: | ---: | ---: |
| 1 | `cell_mean_speed` | 0.1896 | 0.1882 | -0.7% |
| 1 | `cell_straightness` | 0.3225 | 0.3183 | -1.3% |
| 1 | `cell_angular_velocity_rms` | 1.3013 | 1.3161 | +1.1% |
| 2 | `cell_mean_speed` | 0.3767 | 0.3695 | -1.9% |
| 2 | `cell_straightness` | 0.1775 | 0.1673 | -5.7% |
| 2 | `cell_angular_velocity_rms` | 2.3178 | 2.3753 | +2.5% |
| 3 | `cell_mean_speed` | 0.4151 | 0.4242 | +2.2% |
| 3 | `cell_straightness` | 0.1559 | 0.1309 | -16.0% |
| 3 | `cell_angular_velocity_rms` | 5.7270 | 7.1107 | +24.2% |

`n=1,2` は主要特徴の変化が比較的小さい。`n=3` は平均速度を維持する一方，直進性低下と角速度増加が大きいため，training dataとしての自然さをvisual reviewで重点確認した。

## Visual Review Result

2026-07-16 にユーザーが次のreplayを確認し，`n=1,2,3` の定性判断を承認した。

```bash
uv run python scripts/01_simulate_swimming/render_shape_stability_grid_replay.py config=conf/phase2_multi_run/flagella_count_behavior_v1.yaml overwrite=true
```

確認先:

- `outputs/phase2_multi_run/flagella_count_behavior_v1/replay/grid_swim3d_page01.mp4` から `grid_swim3d_page04.mp4`
- 各 `grid_swim3d_page*_final.png`
- `shape_stability_metrics.csv` と `shape_stability_metrics.png`

確認結果:

- `n=1,2,3`: body deformation，helical shape preservation，swimming-like propulsionにtraining candidateから外す問題は見られなかった。
- `n=3`: v0比の直進性低下と角速度増加を踏まえて重点確認し，定性判断は問題なしとした。
- `n=4`: 大きな形状崩壊はなく全体形状は良好だが，べん毛回転が滑らかな等速回転に見えない。strict QC failureも3/9に残るためdiagnostic-onlyを維持する。

## n=4 Rotation Smoothness

v1は `motor.enable_switching=false`，Brownian motion無効，torque ramp無効であり，各べん毛へ与える `motor_torque_Nm` は全stepで `2.0e-20 N m` だった。したがって，`n=4` の非等速回転は入力torque値の時間変動ではない。

`flag_helix_spin_phase_deg`を0.1 s以降に20 ms区間で差分し，回転速度の絶対値に対する変動係数と回転方向整合率を比較した。

| n_flagella | mean rate variation coefficient | mean direction consistency | max flag bond rel err | max torsion err |
| ---: | ---: | ---: | ---: | ---: |
| 1 | 0.219 | 1.000 | 0.174 | 12.3 deg |
| 2 | 0.595 | 0.970 | 0.409 | 17.6 deg |
| 3 | 1.293 | 0.927 | 0.853 | 34.8 deg |
| 4 | 1.746 | 0.848 | 1.151 | 157.8 deg |

入力torque総量は各べん毛で一定でも，`diffusive` profileによる作用位置，RPY mobility，spring / bend / torsion，body反作用は時々刻々の形状に依存する。べん毛数の増加に伴ってこれらのelastic / hydrodynamic load変動と局所変形が増え，回転速度変動が大きくなったと推測する。`n=4` の `flag_helix_spin_fit_r2` は平均0.937，最小0.401まで低下しているため，一部には螺旋変形による位相推定不安定も含まれる。寄与の切り分けと改善はIssue #124で扱う。

## Interpretation

full v1 raw datasetとanalysis datasetは生成済みであり，自動QCとuser visual reviewの両方に基づいてtraining candidateを `n=1,2,3` に確定した。`n=4` はdiagnostic-only，`n=5,6` は対象外とする。Issue #119の完了条件を満たしたため，Phase3/4 training dataset作成へhandoffする。
