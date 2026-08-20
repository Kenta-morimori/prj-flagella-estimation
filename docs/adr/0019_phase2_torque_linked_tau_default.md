# ADR 0019: Phase 2 torque-linked τ を既定とする

* Status: Accepted
* Date: 2026-08-20

## Context

2010 project profileだけが暗黙に固定 `tau_s=1 s` を使っており、同じ物理トルクを変えても
`dt_star` と実時間の対応が profile により異なっていた。

## Decision

`time.scale_policy` 未指定時は、2010/2015、paper/project を問わず
`reference_torque` を使う。標準設定では `motor.torque_Nm` と
`motor.reference_torque_Nm` は同値とし、force-scale override は禁止する。

過去結果を再現する config と固定τ比較 control は
`legacy_fixed_tau_s_1` を明示する。torque/reference torque の不一致または force override は
`motor.allow_reference_torque_mismatch=true` を持つ比較実験だけで許可する。

## Consequences

新規 run の manifest は torque、reference torque、material force scale、`tau_s`、内部刻みを
同時に記録する。τ policy変更だけでは dataset version や training candidate を採択しない。
