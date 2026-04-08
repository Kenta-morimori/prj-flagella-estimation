# PR #42: Phase2/Phase3 baseline2 vs after_split_final 比較（2026-04-08）

## 比較条件
- baseline: outputs/phase_diagnostics_baseline2
- after: outputs/phase_diagnostics_after_split_final
- 対象 phase/time: Phase2・Phase3 × 0.01 s, 0.1 s, 2.0 s
- 判定基準: finite 完走のみでは完了扱いにせず、baseline2 からの局所誤差改善を主軸に判定

## 1ページ比較表（baseline -> after）

| phase | t_s | pos_all_finite | local_attach_first_rel_err | local_first_second_rel_err | flag_bond_rel_err_max | local_first_torsion_err_deg | motor_attach_force_norm | motor_first_force_norm | motor_second_force_norm |
|---|---:|---|---|---|---|---|---|---|---|
| phase2 | 0.010 | True -> True | 7200.505% -> 1442.113% | 4261.267% -> 1181.928% | 4261.267% -> 1181.928% | NaN -> NaN | 6.362e-16 -> 1.992e-15 | 2.302e-16 -> 1.141e-15 | 4.060e-16 -> 8.510e-16 |
| phase2 | 0.100 | True -> True | 8746.052% -> 2642.267% | 4124.919% -> 1156.687% | 4124.919% -> 1171.926% | NaN -> NaN | 7.310e-16 -> 2.364e-15 | 9.201e-17 -> 1.870e-16 | 6.390e-16 -> 2.177e-15 |
| phase2 | 1.999 | True -> True | 5067.795% -> 2085.042% | 863.768% -> 907.639% | 5337.556% -> 1602.929% | NaN -> NaN | 4.097e-16 -> 8.346e-16 | 1.115e-16 -> 1.183e-15 | 5.213e-16 -> 3.488e-16 |
| phase3 | 0.010 | True -> True | 11784.291% -> 11830.985% | 6305.427% -> 6337.213% | 17602.812% -> 17605.631% | 25.138 -> 25.129 | 2.353e-16 -> 2.350e-16 | 8.278e-17 -> 8.285e-17 | 1.525e-16 -> 1.521e-16 |
| phase3 | 0.100 | True -> True | 15777.011% -> 15999.423% | 5187.382% -> 5080.634% | 15614.998% -> 15601.708% | 20.580 -> 17.810 | 1.931e-16 -> 2.020e-16 | 1.501e-17 -> 2.812e-17 | 2.081e-16 -> 2.301e-16 |
| phase3 | 1.999 | True -> True | 29244.593% -> 30719.632% | 2603.342% -> 2667.827% | 9916.905% -> 9781.327% | 9.068 -> 10.580 | 1.281e-16 -> 1.232e-16 | 4.087e-16 -> 4.124e-16 | 5.369e-16 -> 5.356e-16 |

## DoD 判定
- Phase B (Phase2): **達成**
  - 0.1 s で主要局所誤差が baseline2 比で明確に低下
  - attach-first: 8746.052% -> 2642.267%
  - first-second: 4124.919% -> 1156.687%
  - bond max: 4124.919% -> 1171.926%
  - 全対象時刻で pos_all_finite=True を維持

- Phase C (Phase3): **未達**
  - 0.1 s で first-second / bond は改善する一方、attach-first は悪化
  - attach-first: 15777.011% -> 15999.423% (悪化)
  - first-second: 5187.382% -> 5080.634% (改善)
  - bond max: 15614.998% -> 15601.708% (微改善)
  - 2.0 s では attach-first / first-second / torsion が再悪化しており、局所誤差改善の一貫性が不足

## 判定根拠（短文）
- finite 完走は両 phase で達成しているが、DoD はそれ単独では満たさない。
- Phase B は 0.1 s の主要局所誤差3指標で baseline2 比の有意低下を確認できるため DoD 達成。
- Phase C は attach-first の悪化が残り、時刻横断で改善が安定しないため DoD 未達。

## 次アクション（1つ）
- **compute_motor_forces() の attach-first 側 preload 成分のみを再設計**する。
  - 方針: first-second 側は現行維持、attach-first 側の軸方向 preload 係数と投影方向だけを探索して、0.1 s の attach-first を優先改善する。
  - 制約: torque はこれ以上下げない。local constraint は固定。主戦場は compute_motor_forces()。
