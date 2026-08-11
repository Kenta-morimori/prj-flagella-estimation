# Phase 2 Issue #183: reference torque 比較契約

この契約は #61 と #184 の入力であり、torque・`dt_star`・dataset v2・2015 profile supported statusを採択しない。正式な比較policyは [ADR 0016](../adr/0016_phase2_reference_torque_comparison_policy.md) を正本とする。

`motor.torque_Nm` は各flagellar motorの駆動、`motor.reference_torque_Nm` は時間・物性の参照、`torque_for_forces_Nm` は実際の物性係数scaleである。force overrideがなければ後者はreference torqueの絶対値になる。manifestの `time.material_coefficients` は spring `H`、bend/torsion/hook stiffness、repulsion `A` の実効基準係数を保存する。

`fixed-reference` と `tracking-reference`、`same-real-time` と `same-dimensionless-time` は必ず別 condition とする。前者は同一物性での駆動感度、後者は時間と物性を連動させた相似候補である。同一無次元時間を、同一実時間でのstep数・wall time・計算効率の比較に使わない。

2010 projectは `tau_s=1 s` 固定であり、tracking-referenceは時間相似でない。2015 projectは `tau_s=eta*b^3/abs(reference_torque_Nm)` である。#61 の `dt_star` 評価は `fixed-reference` / `same-real-time` に限定して、seed・geometry・物性・per-flag torque policyを固定する。

## 実行準備（simulationは起動しない）

```bash
uv run python scripts/01_simulate_swimming/plan_reference_torque_comparison.py \
  --config conf/phase2_reference_torque/2015_project.yaml \
  --output /private/tmp/issue183-2015-plan.json
```

生成planの各conditionは、effective torque、time manifest、物性係数、実行commandを含む。0.5秒run、dataset採択、2015 supported昇格はこのIssueの範囲外である。
