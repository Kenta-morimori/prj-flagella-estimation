# ADR 0003: Phase 2 local twist torque transmission

- status: proposed
- date: 2026-06-02
- scope: Phase 2.6 P2-6-007

## Context

Phase 2.6では、`triplet` motor が root 近傍に torque を入れても、螺旋全体の安定した一方向回転へ十分に伝わらないことを確認した。

この問題は、時間積分を順次進めていないからではない。bead positions は毎step更新されるため、位置変化は周辺beadへ影響する。しかし、現行モデルには「segment が自分の軸まわりに何度回っているか」という状態量がない。そのため、rootで入ったねじりを、隣接segmentへねじれとして渡す内部経路がない。

現行の torsion force は、4つのbead位置からdihedral angleを計算し、形状が崩れたら戻す力である。これは螺旋形状維持には効くが、root motor torque を segment ごとの軸まわり回転として保存・伝搬する仕組みではない。

## Decision

次の物理モデル候補として、local twist state を導入する。

考え方は以下。

1. flagellum の各 segment に、軸まわりの向き `orientation_i` を持たせる。
2. 隣接segment同士の向きの差を `local_twist_i` として計算する。
3. `local_twist_i` が自然な値 `rest_twist_i` からずれた場合、twist potential で戻そうとする。
4. root motor torque は root 側のorientationへ入力する。
5. 隣接segment間のtwist potentialにより、root側のねじれが先端側へ順に伝わる。

単純化した式は以下。

```text
local_twist_i = orientation_{i+1} - orientation_i - rest_twist_i
U_twist_i = 0.5 * k_twist * local_twist_i^2
torque_i = -dU_twist / d orientation_i
```

これは、べん毛全体に1つのねじれ量を与える方法ではない。連結するsegment同士の相対的な向きから、局所的なねじれを計算する方法である。

## Consequences

- `triplet` motor の root torque を、bead位置の局所変形ではなく、segment twist state として受け取れる。
- root で入ったねじれが、twist potential により隣接segmentへ伝わる。
- `axial_torque_flux_probe` のように離れたbeadへ直接接線力を分布させる必要を減らせる可能性がある。
- 一方で、現行の bead-position-only model からは明確な拡張になる。
- orientation state の時間積分、twist damping、bead positions への力/torque変換を追加する必要がある。

## Initial implementation plan

実装は段階化する。

1. **診断状態の追加**
   - flagellum segment ごとに `orientation_i` と `local_twist_i` を保持・出力できるようにする。
   - この段階では既存の力学挙動を変えず、状態量の定義とログを固定する。

2. **root torque 入力の追加**
   - `triplet` motor torque を root orientation へ入力する。
   - local twist がroot側から増えることを確認する。

3. **twist relaxation の追加**
   - 隣接segment間のtwist potentialにより、local twist が先端側へ伝わることを確認する。
   - まずはorientation state上の伝搬を確認し、bead positionsへの作用は分けて検証する。

4. **bead force / torque への変換**
   - local twist torque を、segment周辺beadsへの接線力または等価force coupleとして変換する。
   - `axial_torque_flux_probe` の0.5 s pass条件を比較基準にする。

## Acceptance criteria

- `time.dt_star=1.0e-4`
- `duration_s >= 0.5`
- `motor.force_distribution=triplet` または root motor由来のlocal twist mode
- `net_abs_flag_helix_spin_revolutions >= 1.0`
- `flag_helix_spin_direction_consistency >= 0.5`
- `hook_len_rel_err_max <= 0.5`
- `flag_bond_rel_err_max <= 0.25`
- `flag_bend_err_max_deg <= 30`
- `flag_torsion_err_max_deg <= 60`

## Relation to probes

`distributed_flagellum` は、螺旋全体へtorqueが届いた場合の上限診断である。

`axial_torque_flux_probe` は、root起点のtorque fluxをbead-position-only modelで近似する短期的な妥協案である。

local twist transmission は、これらのprobeより物理的に根本側へ戻した本命案である。まずprobeと同等の回転・形状維持を達成できるかで評価する。
