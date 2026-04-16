# PR #42 継続: Phase C attach-first preload 再設計（2026-04-08）

## 変更点
- 対象: `compute_motor_forces()` の attach-first preload 係数のみ
- 変更: `preload_attach_first_scale: -0.01 -> +0.01`
- 固定: canonical torque=4e-21, local constraint固定, projection再導入なし

## 比較対象
- baseline2: `outputs/phase_diagnostics_baseline2`
- previous after: `outputs/phase_diagnostics_after_split_final`
- current after: `outputs/phase_diagnostics_after_attach_plus001_d001|d01|d20`

## Phase C 比較表（baseline2 / previous after / current after）

| t_s | pos_all_finite | local_attach_first_rel_err | local_first_second_rel_err | flag_bond_rel_err_max | local_first_torsion_err_deg | motor_attach_force_norm | motor_first_force_norm | motor_second_force_norm |
|---:|---|---|---|---|---|---|---|---|
| 0.010 | True / True / True | 11784.291% / 11830.985% / 11712.290% | 6305.427% / 6337.213% / 6317.834% | 17602.812% / 17605.631% / 17637.179% | 25.138 / 25.129 / 26.344 | 2.353e-16 / 2.350e-16 / 2.435e-16 | 8.278e-17 / 8.285e-17 / 8.572e-17 | 1.525e-16 / 1.521e-16 / 1.578e-16 |
| 0.100 | True / True / True | 15777.011% / 15999.423% / 15712.137% | 5187.382% / 5080.634% / 5256.111% | 15614.998% / 15601.708% / 15607.222% | 20.580 / 17.810 / 15.192 | 1.931e-16 / 2.020e-16 / 1.789e-16 | 1.501e-17 / 2.812e-17 / 1.546e-17 | 2.081e-16 / 2.301e-16 / 1.942e-16 |
| 2.000 | True / True / True | 29244.593% / 30719.632% / 30485.515% | 2603.342% / 2667.827% / 2637.548% | 9916.905% / 9781.327% / 9894.331% | 9.068 / 10.580 / 10.409 | 1.281e-16 / 1.232e-16 / 1.277e-16 | 4.087e-16 / 4.124e-16 / 4.403e-16 | 5.369e-16 / 5.356e-16 / 5.679e-16 |

## 判定（0.1 s hard gate）
- local_attach_first_rel_err: 15777.011% -> 15999.423% -> 15712.137%
  - previous after 比で改善、baseline2 比でも僅かに改善
- local_first_second_rel_err: 5187.382% -> 5080.634% -> 5256.111%
  - previous after 比で悪化
- flag_bond_rel_err_max: 15614.998% -> 15601.708% -> 15607.222%
  - previous after 比で微悪化
- local_first_torsion_err_deg: 20.580 -> 17.810 -> 15.192
  - 改善
- 結論: attach-first 主指標は改善したが、補助指標（first-second/bond）の同時改善が不十分で、Phase C DoD は未達。

## 次アクション（1つ）
- `compute_motor_forces()` の attach-first preload を定数値から局所状態依存へ移行し、0.1 s 時点で attach-first を維持しつつ first-second 悪化を抑える重み関数を導入する（torque/local constraint は固定のまま）。
