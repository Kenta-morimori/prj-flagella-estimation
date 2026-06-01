# Phase 2.6 triplet torque transmission diagnostics

## 目的

`distributed_flagellum` で綺麗に回転することを確認した後、同じブランチで `triplet + hook spring fix` を再評価し、root 方位の回転が螺旋全体へどれだけ伝わるかを定量化する。

## 実施内容

- `summarize_single_flagellum_helix_retention()` に root-to-helix 伝達診断を追加した。
  - `net_abs_flag_root_revolutions`
  - `signed_flag_root_rate_hz`
  - `helix_to_root_net_rotation_ratio`
- `dt_star=1.0e-4`, 0.25 s 条件で `triplet + hook spring fix` を再評価した。
- `docs/phase2/phase2_6_torque_transmission_notes.md` に、ねじれ・回転自由度不足の整理と次タスクを記録した。
- `docs/phase2/phase2_tasks.md` に P2-6-007 を追加した。
- ユーザー目視レビュー完了を受け、`20260601_134936_phase2_p2-6_distributed_spin/review_result.json` を PASS に更新した。

## triplet 再評価結果

共通条件:

- `stub_mode=full_flagella`
- `n_flagella=1`
- `dt_star=1.0e-4`
- `duration_s=0.25`
- `motor.force_distribution=triplet`
- hook spring force 復元後

| torque | local scales `(hook,spring,bend,torsion)` | first fail | root net rev | helix net rev | transfer ratio |
| ---: | --- | --- | ---: | ---: | ---: |
| `2.0e-20` | `(1,1.2,1,1)` | `motor_no_rotation` | `0.0665` | `0.00340` | `0.0510` |
| `2.0e-20` | `(1,1,1,1)` | `motor_no_rotation` | `0.0665` | `0.00340` | `0.0510` |
| `2.5e-20` | `(1,1.2,1,1)` | `motor_no_rotation` | `0.0465` | `0.0448` | `0.962` |
| `2.5e-20` | `(4,2,2,2)` | `motor_no_rotation` | `0.8796` | `0.00139` | `0.00158` |

## 解釈

`triplet` は root 近傍へ torque を入れているが、現行 bead-position-only model には material frame / segment twist / axial torque flux がない。そのため root 方位の運動が螺旋全体の一方向 net 回転へ安定して伝わらない。

`distributed_flagellum` は diagnostic upper-bound として残し、次は root torque を chain へ伝える物理自由度の設計を検討する。
