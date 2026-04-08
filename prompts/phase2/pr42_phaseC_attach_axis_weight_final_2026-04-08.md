# PR #42: Phase C attach preload final comparison（2026-04-08）

## 変更点
- attach-first preload を定数値から body axis 依存の重み関数へ変更
- body axis 依存量: `motor_axis_vs_rear_direction_angle_deg` 相当
- first-second 側 preload は `+0.01` のまま維持
- canonical torque = `4e-21`、local constraint は固定

## 比較条件
- baseline: `outputs/phase_diagnostics_baseline2`
- previous after: `outputs/phase_diagnostics_after_split_final`
- current after: `outputs/phase_diag_attach_axis_weight_final`
- 対象 phase/time: Phase2・Phase3 × 0.01 s, 0.1 s, 2.0 s

## 1ページ比較表（baseline / previous after / current after）

| phase | t_s | pos_all_finite | local_attach_first_rel_err | local_first_second_rel_err | flag_bond_rel_err_max | local_first_torsion_err_deg | motor_attach_force_norm | motor_first_force_norm | motor_second_force_norm |
|---|---:|---|---|---|---|---|---|---|---|
| phase2 | 0.010 | True / True / True | 7200.505% / 1442.113% / 1451.871% | 4261.267% / 1181.928% / 1186.232% | 4261.267% / 1181.928% / 1186.232% | NaN / NaN / NaN | 6.362e-16 / 1.992e-15 / 1.992e-15 | 2.302e-16 / 1.141e-15 / 1.140e-15 | 4.060e-16 / 8.510e-16 / 8.517e-16 |
| phase2 | 0.100 | True / True / True | 8746.052% / 2642.267% / 2615.712% | 4124.919% / 1156.687% / 1096.317% | 4124.919% / 1171.926% / 1225.591% | NaN / NaN / NaN | 7.310e-16 / 2.364e-15 / 4.932e-15 | 9.201e-17 / 1.870e-16 / 8.345e-17 | 6.390e-16 / 2.177e-15 / 4.849e-15 |
| phase2 | 2.000 | True / True / True | 5067.795% / 2085.042% / 3805.839% | 863.768% / 907.639% / 910.741% | 5337.556% / 1602.929% / 1611.603% | NaN / NaN / NaN | 4.097e-16 / 8.346e-16 / 4.150e-16 | 1.115e-16 / 1.183e-15 / 1.869e-16 | 5.213e-16 / 3.488e-16 / 2.283e-16 |
| phase3 | 0.010 | True / True / True | 11784.291% / 11830.985% / 11830.932% | 6305.427% / 6337.213% / 6337.188% | 17602.812% / 17605.631% / 17605.629% | 25.138 / 25.129 / 25.129 | 2.353e-16 / 2.350e-16 / 2.350e-16 | 8.278e-17 / 8.285e-17 / 8.285e-17 | 1.525e-16 / 1.521e-16 / 1.521e-16 |
| phase3 | 0.100 | True / True / True | 15777.011% / 15999.423% / 15890.223% | 5187.382% / 5080.634% / 5041.721% | 15614.998% / 15601.708% / 15595.971% | 20.580 / 17.810 / 18.521 | 1.931e-16 / 2.020e-16 / 1.856e-16 | 1.501e-17 / 2.812e-17 / 2.030e-17 | 2.081e-16 / 2.301e-16 / 2.057e-16 |
| phase3 | 2.000 | True / True / True | 29244.593% / 30719.632% / 29452.406% | 2603.342% / 2667.827% / 2664.555% | 9916.905% / 9781.327% / 9787.495% | 9.068 / 10.580 / 10.573 | 1.281e-16 / 1.232e-16 / 1.251e-16 | 4.087e-16 / 4.124e-16 / 4.029e-16 | 5.369e-16 / 5.356e-16 / 5.280e-16 |

## Phase C 判定
- `local_attach_first_rel_err`: baseline2 よりも改善、previous after よりも改善
- `local_first_second_rel_err`: baseline2 よりは悪化、previous after よりも悪化
- `flag_bond_rel_err_max`: baseline2 よりは悪化、previous after よりも悪化
- `local_first_torsion_err_deg`: previous after より改善
- 結論: Phase C は still 未達。attach-first の改善は維持できたが、first-second / bond の悪化を抑え切れていない。

## DoD までの距離
- 0.1 s での attach-first: 15777.011% -> 15999.423% -> 15890.223%
- 0.1 s での first-second: 5187.382% -> 5080.634% -> 5041.721%
- 0.1 s での bond max: 15614.998% -> 15601.708% -> 15595.971%

## 次アクション（1つ）
- `compute_motor_forces()` の attach-first 重み関数に対して、phase3 側のみさらに局所的に効く第2の1変数条件を追加するのではなく、まずは現行の body axis 依存関数の係数を 1 回だけ微調整する。
- その後も first-second / bond が previous after を上回れないなら、attach-first preload の影響範囲を縮める方向で検討する。
