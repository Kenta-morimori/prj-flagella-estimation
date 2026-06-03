# ADR 0003: Phase 2 local twist transmission probe history

- status: superseded
- date: 2026-06-02
- scope: Phase 2.6 P2-6-007
- superseded_by:
  - `docs/adr/0002_phase2_torque_transmission_probes.md`
  - `docs/adr/0004_phase2_material_frame_twist_transmission.md`

## Context

このADRは当初、`local_twist_transmission_probe` の設計と material frame / segment twist による完全物理実装方針を同時に記録していた。

ユーザー確認により、本PRの scope は「probe成功まで」とし、完全物理実装は次タスクへ分離する方針になった。そのため、probe の実装検証結果は ADR 0002 に統合し、完全物理実装の方針は ADR 0004 に分離する。

## Historical decision

P2-6-007では、実験用 `motor.force_distribution=local_twist_transmission_probe` を追加した。

この mode は、flagellum segment ごとに内部 `orientation_i` を持ち、root motor torque により root側の orientation activity を進める。その activity を diffusion/relaxation で先端側へ伝え、bead接線forceの重みとして使う。

代表条件では以下を確認した。

- `motor.force_distribution=local_twist_transmission_probe`
- `motor.torque_Nm=2.0e-20`
- `time.dt_star=1.0e-4`
- `duration_s=0.5`
- `local_spring_scale=1.2`
- `net_abs_flag_helix_spin_revolutions=1.04745`
- `flag_helix_spin_direction_consistency=0.81507`
- `local_twist_tip_activity_ratio=0.34836`

## Supersession

この結果は、root側orientation activityが先端側へ伝われば、単一べん毛で螺旋形状を保ちつつnet回転できることを示す。ただし、local twist potential から局所的な bead force / force couple を導出していないため、完全な material frame / segment twist 物理モデルではない。

今後は以下のように扱う。

- probe の比較・結果: ADR 0002
- 完全物理実装の設計方針: ADR 0004
