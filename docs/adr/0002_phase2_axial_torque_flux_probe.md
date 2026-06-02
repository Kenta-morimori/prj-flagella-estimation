# ADR 0002: Phase 2 axial torque flux probe

- status: accepted
- date: 2026-06-02
- scope: Phase 2.6 P2-6-007

## Context

Phase 2.6 では、`triplet + hook spring fix` が root 方位の net 回転を作っても、螺旋全体の net 回転へ十分に伝わらないことを確認した。

代表的な `triplet` 条件では、`net_abs_flag_root_revolutions=0.8796` に対して `net_abs_flag_helix_spin_revolutions=0.00139` であり、`helix_to_root_net_rotation_ratio=0.00158` だった。

この失敗は、hook spring force だけの問題ではない。現行 `triplet` は `attach-first-second` の局所3点へ力カップルを入れるが、現行 bead-position-only model には segment の軸まわり姿勢や、root から先端側へ torque flux を輸送する状態量がない。そのため root 方位 proxy は動いても、螺旋全体の一方向接線速度へ変換されにくい。

一方で、`distributed_flagellum` は螺旋全体へ直接 torque を分布させるため、単一べん毛の螺旋形状維持と net 回転を両立できた。ただし、これは root motor torque の伝搬モデルではなく、diagnostic upper-bound である。

## Decision

`motor.force_distribution=axial_torque_flux_probe` を追加する。

この mode は、root bead を原点、flagellum の root-to-tip 主軸を torque flux 軸とし、root から先端へ緩やかに減衰する接線力を flagellum beads に分布させる。同じ軸まわりの反作用 torque は body 側へ与える。

この mode は material frame / segment twist を明示する本実装ではない。現行 bead-position-only model で「軸方向 torque flux 近似が有効か」を確認するための probe として扱う。

`distributed_flagellum` との違いは、先端側の減衰だけではない。`distributed_flagellum` は flagellum 全体の重心まわりに、全体へほぼ一様な接線driveを入れる。一方、`axial_torque_flux_probe` は root beadを原点にし、root側を強く先端側を弱くする。どちらも離れたbeadへ直接力を入れる近似だが、`distributed_flagellum` より root motor torque の伝搬解釈に近い比較用probeとして扱う。

## Consequences

- `triplet` 既定値は変更しない。
- `distributed_flagellum` は診断用 upper-bound として残す。
- `axial_torque_flux_probe` は、`distributed_flagellum` より root motor 由来に近い probe だが、離れた bead に直接力を入れる近似である。
- 参照論文モデルからの拡張であり、最終的な物理モデルとして採用する前に material frame / segment twist との比較が必要である。
- 本PRでは、短期比較用probeを `axial_torque_flux_probe`、上限診断を `distributed_flagellum` として分ける。P2-6-007後半では、内部orientation stateを持つ `local_twist_transmission_probe` を妥協案候補として優先する判断に更新した。詳細は `docs/adr/0003_phase2_local_twist_transmission.md` に記録する。

## Verification

短時間条件:

- `motor.force_distribution=axial_torque_flux_probe`
- `torque=2.0e-20 N m`
- `time.dt_star=1.0e-4`
- `duration_s=0.05`
- `local_hook_scale=1.0`
- `local_spring_scale=1.2`
- `local_bend_scale=1.0`
- `local_torsion_scale=1.0`

この条件で `net_abs_flag_helix_spin_revolutions=0.11566`, `flag_helix_spin_direction_consistency=1.0`, `helix_to_root_net_rotation_ratio=2.81489` を確認した。

長時間条件:

- `motor.force_distribution=axial_torque_flux_probe`
- `torque=2.0e-20 N m`
- `time.dt_star=1.0e-4`
- `duration_s=0.5`
- `local_hook_scale=1.0`
- `local_spring_scale=1.2`
- `local_bend_scale=1.0`
- `local_torsion_scale=1.0`

この条件で `net_abs_flag_helix_spin_revolutions=1.08546`, `flag_helix_spin_direction_consistency=1.0`, `hook_len_rel_err_max=0.39986`, `flag_bond_rel_err_max=0.17206`, `flag_bend_err_max_deg=1.12490`, `flag_torsion_err_max_deg=3.97803` を確認した。
