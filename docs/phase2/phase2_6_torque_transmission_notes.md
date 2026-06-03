# Phase 2.6 torque transmission notes

## 背景

`distributed_flagellum` により、単一 full flagellum は `dt_star=1.0e-4`, 0.5 s 条件で螺旋形状を保ちつつ net 1回転以上できた。これは「トルクが螺旋全体へ伝われば綺麗に回る」ことを示す診断結果である。

一方で、従来の `triplet` motor は `attach-first-second` の3点に局所的な力カップルを与える。hook spring force を復元しても、root 方位の運動が螺旋全体の net 回転へ十分に伝わらない条件が残る。

## 再評価結果

共通条件:

- `stub_mode=full_flagella`
- `n_flagella=1`
- `dt_star=1.0e-4`
- `duration_s=0.25`
- `motor.force_distribution=triplet`
- hook spring force 復元後

| torque | local scales `(hook,spring,bend,torsion)` | first fail | root net rev | helix net rev | transfer ratio | 備考 |
| ---: | --- | --- | ---: | ---: | ---: | --- |
| `2.0e-20` | `(1,1.2,1,1)` | `motor_no_rotation` | `0.0665` | `0.00340` | `0.0510` | 形状は保つが回転が弱い |
| `2.0e-20` | `(1,1,1,1)` | `motor_no_rotation` | `0.0665` | `0.00340` | `0.0510` | spring scale の影響は小さい |
| `2.5e-20` | `(1,1.2,1,1)` | `motor_no_rotation` | `0.0465` | `0.0448` | `0.962` | helix phase は往復揺れが大きく direction consistency `0.055` |
| `2.5e-20` | `(4,2,2,2)` | `motor_no_rotation` | `0.8796` | `0.00139` | `0.00158` | root 方位は動くが螺旋全体へほぼ伝わらない |

## 解釈

現行 `triplet` motor は root 周辺に torque を入れるが、現行の bead-position-only model には segment ごとの material frame や明示的な twist degree of freedom がない。そのため、root 方位の運動が局所変形、弾性力との釣り合い、またはフィット jitter として消え、螺旋全体の一方向 net 回転へ伝わりにくい。

`distributed_flagellum` は、この不足を補う diagnostic upper-bound である。これは root motor の厳密モデルではなく、「もし螺旋全体へ torque が伝わればどう見えるか」を示す比較対象として扱う。

## 次の設計候補

1. **material frame / segment twist の導入**
   - 各 segment に断面の向きまたは material frame を持たせる。
   - root torque を twist として chain に伝える。
   - 物理的には最も明確だが、状態変数と積分器の拡張が必要。

2. **segment torsional torque flux**
   - 既存 bead positions に加え、隣接 segment 間へ軸方向 torque flux を近似的に分配する。
   - `distributed_flagellum` より root 由来に近いが、実装と安定性の検証が必要。

3. **quasi-rigid helical body approximation**
   - 単一 flagellum の螺旋形状をほぼ剛体として扱い、root torque で全体を回す。
   - ML 用データ生成には実用的だが、bead-spring 論文モデルからは離れる。

## 評価基準

- `time.dt_star=1.0e-4`
- `duration_s >= 0.5`
- `net_abs_flag_helix_spin_revolutions >= 1.0`
- `flag_helix_spin_direction_consistency >= 0.5`
- `helix_to_root_net_rotation_ratio` が diagnostic baseline より改善すること
- `hook_len_rel_err_max <= 0.5`
- `flag_bond_rel_err_max <= 0.25`
- `flag_bend_err_max_deg <= 30`
- `flag_torsion_err_max_deg <= 60`
