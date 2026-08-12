# Phase 2 Issue #183: reference torque 比較契約

この契約は #61 と #184 の入力であり、torque・`dt_star`・dataset v2・2015 profile supported statusを採択しない。正式な比較policyは [ADR 0016](../adr/0016_phase2_reference_torque_comparison_policy.md) を正本とする。

`motor.torque_Nm` は各flagellar motorの駆動、`motor.reference_torque_Nm` は時間・物性の参照、`torque_for_forces_Nm` は実際の物性係数scaleである。force overrideがなければ後者はreference torqueの絶対値になる。manifestの `time.material_coefficients` は spring `H`、bend/torsion/hook stiffness、repulsion `A` の実効基準係数を保存する。

時間のuser-facing表記は論文に合わせる。`t` は実時間[s]、`τ` は時間尺度[s]、`Δt` は実時間の積分刻み[s]であり、次が成り立つ。

```text
t / τ = t_star
Δt / τ = dt_star
Δt = dt_star * τ
```

manifestの `time.paper_notation` はこの表記を正本の表示層として保存する。JSON keyはASCIIで、`t`、`tau`、`delta_t`、`t_over_tau`、`delta_t_over_tau` とする。既存の `tau_s`、`duration_tau`、`dt_star`、`dt_internal_s` は、既存runと解析コードとの互換のため残す。

`fixed-reference` と `tracking-reference`、`same-real-time` と `same-dimensionless-time` は必ず別 condition とする。前者は同一物性での駆動感度、後者は時間と物性を連動させた相似候補である。同一無次元時間を、同一実時間でのstep数・wall time・計算効率の比較に使わない。

2010 projectは `tau_s=1 s` 固定であり、tracking-referenceは時間相似でない。2015 projectは `tau_s=eta*b^3/abs(reference_torque_Nm)` である。#61 の `dt_star` 評価は `fixed-reference` / `same-real-time` に限定して、seed・geometry・物性・per-flag torque policyを固定する。

## 実行準備（simulationは起動しない）

各 policy はさらに別々に、`same-real-time`（`time.duration.unit=s`）と `same-dimensionless-time`（`unit=tau`）を実行する。前者だけが #61 の計算効率、wall time、steps/s の主比較である。後者は無次元軌道・離散化の補助比較であり、実時間効率の根拠にしてはならない。

2010 project は `legacy_fixed_tau_s_1` のため、tracking-reference でも `tau_s` は変わらない。これは物性連動を検査する比較であり、2015 のような時間相似を主張できない。plan manifest の `tau_tracks_reference_torque=false` を必ず確認する。

## Config / manifest 契約

入力 config は [2010 project](../../conf/phase2_reference_torque/2010_project.yaml)、[2015 project](../../conf/phase2_reference_torque/2015_project.yaml) である。plan CLI は condition ごとに以下を `campaign_plan.json` へ固定する。

- `comparison_policy`, `time_basis`, `motor_torque_scale`
- 有効な `motor_torque_Nm`, `reference_torque_Nm`, `torque_for_forces_Nm`
- `time`（`duration_s`, `duration_tau`, `tau_s`, `dt_internal_s`, `total_steps` を含む）
- 同等性 flags と実行用 `execution_command`

各 simulation output は従来どおり `manifest.json` と `sim/run_summary.json` を必須とし、campaign の集計は `summary.csv` にする。通常確認で巨大な `step_summary.csv` を直接読まない。必要時だけ bounded inspection を使う。

## 評価指標と受入条件

各 condition は #168 と同じ基本 safety（例外なし、finite、body/non-body shape gate、motor action-reaction diagnostics）を満たすこと。さらに #61 は同一実時間の `total_steps`、`dt_internal_s`、wall time、steps/s、final/trajectory の数値差を比較する。#184 は #61 が許容した `dt_star` だけを用い、torque × `n_flagella` の safety、replay、heatmap、dataset v2 の採択判断を行う。

この Issue の受入条件は次である。

1. fixed / tracking と実時間 / 無次元時間を manifest 上で機械的に区別できる。
2. 同一実時間の比較で物性・seed・geometry・`dt_star` 以外を揃え、差分を明記する。
3. #61 / #184 は policy 間の軌道差を「数値誤差」や「同一物理系の改善」と解釈しない。
4. Stage A output を再計算・書換え・採択根拠への昇格をしない。

## scale の関係と hook-first diagnostic

`motor_torque_scale` は基準 `reference_torque_Nm` に掛ける倍率である。例えば2010 projectの基準`2.5e-20 N m`では、`0.5 / 1 / 2`は`1.25e-20 / 2.5e-20 / 5.0e-20 N m`を表す。

物性の全体基準は `torque_for_forces_Nm` で、実装上は次を比例させる。

| quantity | base relation |
| --- | --- |
| spring coefficient `H` | `H_over_T_over_b * torque_for_forces_Nm / b` |
| flag/body bend coefficient | `kb_over_T * torque_for_forces_Nm` |
| torsion coefficient | `kt_over_T * torque_for_forces_Nm` |
| hook bend coefficient | `hook.kb_over_T * torque_for_forces_Nm` |
| segment repulsion amplitude | `A_ss_over_T * torque_for_forces_Nm` |

従って`fixed-reference`はこれらを固定し、`tracking-reference`はreference torqueとともに全体を比例させる。`motor.local_hook_scale`はこの全体scaleとは別で、motor-on時にhook forceだけへ掛かる局所倍率である。`local_spring_scale` / `local_bend_scale` / `local_torsion_scale`も全体ではなく根元近傍の対応する項だけに掛かる。

2010 tracking-reference scale=2ではhookが最初にFAILしたため、次の短時間diagnosticを用意する。[hook diagnostic config](../../conf/phase2_multi_run/reference_torque_2010_tracking_scale2_hook_diagnostic.yaml)は`local_hook_scale=1,1.25,1.5,2`だけを変え、その他を固定する。hook PASSだけでは候補とせず、flag bond/bend/torsionを含むshape gateが全期間PASSした場合だけ、次段階の同一`1 s`確認へ進める。hook scaleでflag failureが現れる、または短時間PASSでも1 sでFAILする場合は、hook単独解決とは扱わない。

### 2010 project の先行実行結果（2026-08-11）

`same-real-time=1 s`、`dt_star=1e-4`、10,000 steps、default seedの5条件をユーザーが実行した。有限値は全条件PASSしたが、non-body shapeは次のとおりであった。

| policy | torque scale | shape result | first observed failure |
| --- | ---: | --- | --- |
| fixed-reference | 0.5 | FAIL | `0.9918 s`、flag先行後にhook |
| fixed-reference | 1.0 | PASS | - |
| fixed-reference | 2.0 | FAIL | `0.2698 s`、flag/hook |
| tracking-reference | 0.5 | PASS | - |
| tracking-reference | 2.0 | FAIL | `0.0001 s`、hook |

`tracking-reference` scale=1はfixed-reference scale=1と完全に同一のため実行しなかった。2010 projectは`tau_s=1 s`固定なので、この結果は物性scale連動の先行screenであり、reference torqueに時間scaleも連動する2015相似条件のruntime検証ではない。body diagnosticsは当該runでunavailableであり、body PASSを主張しない。これらの結果はtorque / policyの採択根拠ではない。

## ユーザー実行

まず plan を作る。これは simulation を起動しない。

```bash
uv run python scripts/01_simulate_swimming/plan_reference_torque_comparison.py \
  --config conf/phase2_reference_torque/2015_project.yaml \
  --output /private/tmp/issue183-2015-plan.json
```

生成planの各conditionは、effective torque、time manifest、物性係数、実行commandを含む。0.5秒run、dataset採択、2015 supported昇格はこのIssueの範囲外である。
