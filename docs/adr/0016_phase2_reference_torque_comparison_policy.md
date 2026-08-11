# ADR 0016: Phase 2 reference torque 比較 policy

* Status: Accepted
* Date: 2026-08-12
* Issues: #183, #61, #184, #157

## Context

2010 / 2015 project modelでは、`motor.torque_Nm`（各flagellar motorへ与える符号付き駆動）と、`motor.reference_torque_Nm`（`tau_s` と無次元物性係数の基準）が別の責務を持つ。`torque_for_forces_Nm` は、明示的な `motor.torque_for_forces_override_Nm > 0` があればそれを、なければ reference torque の絶対値を使う。

従って reference torque を motor torque に追従させると、単に `tau_s = eta*b^3/abs(T_ref)` と `dt_s = dt_star*tau_s` が変わるだけではない。spring `H`、bend/torsion/hook stiffness、spring-spring repulsion `A` も `torque_for_forces_Nm` に比例して変わる。実際の係数は run manifest の `time.material_coefficients` に保存する。

## Decision

torque sweep は次の二つを別 campaign として扱う。

| policy | `motor.torque_Nm` | reference / force scale | 比較目的 |
| --- | --- | --- | --- |
| `fixed-reference` | scaleごとに変更 | 基準 `T_ref,0` に固定。force overrideも `T_ref,0` | 同一物性における駆動 torque 感度 |
| `tracking-reference` | scaleごとに変更 | `abs(motor torque)` に連動。force override は0 | 無次元化に沿った相似候補（時間・物性の同時変更） |

各 policy を `same-real-time` と `same-dimensionless-time` に直交させる。#61 の `dt_star` 比較と性能主張は、同一 seed、geometry、物性、per-flag torque policy、duration **を固定した `fixed-reference` / `same-real-time`** 内でのみ行う。`same-dimensionless-time` は無次元軌道を比べる補助比較であり、実時間計算効率を主張しない。

2010 project は `legacy_fixed_tau_s_1` の互換 policy のため、tracking-referenceでも `tau_s=1 s` のままである。物性連動は検査対象だが、2015のような時間相似を主張しない。

## Consequences

* plan config と `campaign_plan.json` は `comparison_policy`、`time_basis`、3種類の torque、`tau_s`、`dt_star`、`dt_internal_s`、`total_steps`、物性係数を保存する。
* #184 は #61 が選んだ `dt_star` を使用し、per-flag torqueを canonicalに保つ。fixed-total-torqueは診断に限る。
* #157 に渡すのはこの比較 contract であり、dataset v2、0.5秒run、2015 supported status の採択ではない。
* #168 / PR #176 Stage A outputs は reference `1.2e-18 N m` 固定の `fixed-reference` 記録として保持し、再実行・遡及解釈しない。

## Evidence

* ADR 0012（time schema と profile別 `tau_s`）
* `conf/phase2_reference_torque/`
* `tests/test_reference_torque_comparison.py`
