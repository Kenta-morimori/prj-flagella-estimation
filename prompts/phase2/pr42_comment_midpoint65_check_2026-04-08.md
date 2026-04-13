# PR #42 中間レポート: body-axis-angle midpoint 65 試行（2026-04-08）

## 目的
- body-axis-angle 重み関数の midpoint を 60 -> 65 に 1 回だけ動かし、Phase B / Phase C への影響を確認した。
- 判定は `step_summary.csv` の実値のみを正とする。

## 比較条件
- current retained candidate: `outputs/phase_diag_attach_axis_weight_final_v2`（midpoint=60）
- trial candidate: `outputs/phase_diag_attach_axis_weight_mid65`（midpoint=65）
- 比較時刻: 0.1 s
- 固定条件: canonical torque = `4e-21`, local constraint 固定, first-second preload 固定

## 0.1 s 比較表
| phase | metric | midpoint=60 | midpoint=65 |
|---|---|---:|---:|
| Phase2 | `local_attach_first_rel_err` | 2649.919% | 1815.522% |
| Phase2 | `local_first_second_rel_err` | 1097.286% | 1140.906% |
| Phase2 | `flag_bond_rel_err_max` | 1224.663% | 1160.827% |
| Phase3 | `local_attach_first_rel_err` | 15950.680% | 16431.162% |
| Phase3 | `local_first_second_rel_err` | 5163.968% | 5122.239% |
| Phase3 | `flag_bond_rel_err_max` | 15606.292% | 15633.602% |

## 判定
- midpoint=65 は Phase B を改善する一方、Phase C の hard gate を壊すため不採用。
- current retained candidate の midpoint=60 を維持する。
- Phase C は 0.1 s で 3 指標すべて previous after より改善しているため、現候補を継続採用する。

## 次アクション
- midpoint は固定し、以後は Phase C の current candidate を比較基準として扱う。
