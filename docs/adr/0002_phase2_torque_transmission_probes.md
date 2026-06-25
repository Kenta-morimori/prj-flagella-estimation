# ADR 0002: Phase 2 torque transmission probes

- status: accepted
- date: 2026-06-02
- scope: Phase 2.6 P2-6-007

## Naming Update

Issue #88 で `axial_torque_flux_probe` / `local_twist_transmission_probe` はコードから削除した。このADR内の probe 名は、P2-6-007当時の診断履歴として読む。旧 `distributed_flagellum` は現行コードでは `root_torque_axis_projection` の deprecated alias である。

## Context

Phase 2.6 では、`triplet + hook spring fix` が root 方位の net 回転を作っても、螺旋全体の net 回転へ十分に伝わらないことを確認した。

代表的な `triplet` 条件では、`net_abs_flag_root_revolutions=0.8796` に対して `net_abs_flag_helix_spin_revolutions=0.00139` であり、`helix_to_root_net_rotation_ratio=0.00158` だった。

この失敗は、hook spring force だけの問題ではない。現行 `triplet` は `attach-first-second` の局所3点へ力カップルを入れるが、現行 bead-position-only model には segment の軸まわり姿勢や、root から先端側へ torque flux を輸送する状態量がない。そのため root 方位 proxy は動いても、螺旋全体の一方向接線速度へ変換されにくい。

`distributed_flagellum` は螺旋全体へ直接 torque を分布させるため、単一べん毛の螺旋形状維持と net 回転を両立できた。ただし、これは root motor torque の伝搬モデルではなく、diagnostic upper-bound である。

## Decision

P2-6-007では、完全物理実装ではなく、torque transmission の有効性を切り分ける probe として以下を追加する。

- `motor.force_distribution=axial_torque_flux_probe`
- `motor.force_distribution=local_twist_transmission_probe`

`axial_torque_flux_probe` は、root bead を原点、flagellum の root-to-tip 主軸を torque flux 軸とし、root から先端へ緩やかに減衰する接線力を flagellum beads に分布させる。同じ軸まわりの反作用 torque は body 側へ与える。

`local_twist_transmission_probe` は、segment ごとの軸まわり orientation state を持ち、root torque による orientation activity を先端側へ拡散・緩和させる。その activity を segment 重みとして使い、root 軸まわりの接線力を beads へ分布させる。

どちらも material frame / segment twist から局所的な force couple を導出する完全物理実装ではない。特に、bead force への最終変換では離れた flagellum bead へ直接接線力を入れる近似が残る。

## Probe の位置づけ

| mode | 目的 | 論文モデルとの距離 | 本PRでの扱い |
| --- | --- | --- | --- |
| `distributed_flagellum` | torque が螺旋全体へ届いた場合の上限診断 | 大。root motor torque の伝搬ではない | ADR 0001に残す |
| `axial_torque_flux_probe` | root起点の torque flux 近似が有効かを確認する | 中-大。material frameなし、非局所接線力あり | 比較用probe |
| `local_twist_transmission_probe` | 内部 orientation state の先端側伝搬が有効かを確認する | 中。内部状態は持つが force変換は近似 | 本PRの主probe |

`distributed_flagellum` と `axial_torque_flux_probe` の違いは、先端側の減衰だけではない。`distributed_flagellum` は flagellum 全体の重心まわりに、全体へほぼ一様な接線driveを入れる。一方、`axial_torque_flux_probe` は root beadを原点にし、root側を強く先端側を弱くする。どちらも離れたbeadへ直接力を入れる近似だが、`axial_torque_flux_probe` の方が root motor torque の伝搬解釈に近い比較用probeである。

`local_twist_transmission_probe` は、`axial_torque_flux_probe` より root torque 伝搬の内部状態を明示する。ただし、local twist potential から bead force を導出していないため、完全な material frame / segment twist 物理モデルとしては扱わない。

## Verification

短時間比較 (`duration_s=0.05`, `time.dt_star=1.0e-4`, `motor.torque_Nm=2.0e-20`, `local_scale=(1,1.2,1,1)`):

| mode | pass | first fail | root net rev | helix net rev | transfer ratio | direction consistency |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| `triplet` | false | `motor_no_rotation` | `0.03789` | `0.00023` | `0.00618` | `0.01690` |
| `axial_torque_flux_probe` | true | `none` | `0.04109` | `0.11566` | `2.81489` | `1.00000` |
| `distributed_flagellum` | true | `none` | `0.05964` | `0.10966` | `1.83882` | `1.00000` |
| `local_twist_transmission_probe` | true | `none` | `0.02072` | `0.12133` | `5.85492` | `1.00000` |

長時間確認 (`duration_s=0.5`, `time.dt_star=1.0e-4`, `motor.torque_Nm=2.0e-20`, `local_scale=(1,1.2,1,1)`):

| mode | pass | helix net rev | direction consistency | hook err | bond err | bend err deg | torsion err deg |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `axial_torque_flux_probe` | true | `1.08546` | `1.00000` | `0.39986` | `0.17206` | `1.12490` | `3.97803` |
| `local_twist_transmission_probe` | true | `1.04745` | `0.81507` | `0.31502` | `0.17314` | `2.52590` | `7.77964` |

`local_twist_transmission_probe` の代表条件では、0.5 s時点で `local_twist_root_orientation_deg=66.809`, `local_twist_tip_orientation_deg=23.273`, `local_twist_abs_max_deg=8.855`, `local_twist_tip_activity_ratio=0.34836` だった。root側に入力した orientation activity が先端側へ部分的に伝わり、その状態に基づくbead force変換で螺旋net 1回転以上とshape gate PASSを両立できることを確認した。

## Consequences

- `triplet` 既定値は変更しない。
- 本PRは probe 成功までを scope とし、完全物理実装は完了条件に含めない。
- `axial_torque_flux_probe` と `local_twist_transmission_probe` は、Phase 2.6 の原因切り分けと次実装の比較基準として残す。
- 完全物理実装では、material frame / segment twist / 局所 force couple を導入し、非局所 force injection を使わないことを目標にする。
- 完全物理実装の方針は ADR 0004 に分離する。
