# Phase 2.4 hook gate 基準

## 目的

Phase 2.4 では、`body + hook` 相当の `minimal_basal_stub` 条件で、motor-on 時の hook 主導 failure を再現可能にする。

この段階では full flagellum の螺旋維持や多本べん毛の束化は対象外とする。対象は basal 近傍、特に attach-first hook link と first-second link の診断である。

## 標準 sweep コマンド

```bash
uv run python -m scripts.01_simulate_swimming.run_motor_scale_sweep \
  --target local_hook_scale \
  --values 1,8 \
  --torques 1.2e-21,4e-21 \
  --duration 0.02 \
  --body-stiffness-scale 50 \
  --output-dir outputs/phase2_4_hook_gate
```

このコマンドは以下を出力する。

- `outputs/phase2_4_hook_gate/local_hook_scale_sweep_summary.csv`
- 条件ごとの `step_summary.csv`
- 条件ごとの `body_constraint_diagnostics.csv`

## Gate 指標

第一の正本は `step_summary.csv` である。

| 指標 | 役割 |
| --- | --- |
| `shape_pass_nonbody` | hook / flagella 側 shape pass |
| `first_fail_category_nonbody` | non-body 側の first-fail category |
| `local_attach_first_rel_err` | attach-first hook link の相対誤差 |
| `hook_len_rel_err_max` | hook 長の最大相対誤差 |
| `hook_angle_err_max_deg` | hook bend 角誤差 |
| `local_first_second_rel_err` | basal first-second link の相対誤差 |

`body_constraint_diagnostics.csv` も記録する。現行の診断 baseline では、高トルク fail 条件で最終 step 時点の body 形状も変形し得る。そのため hook/body の切り分けは以下の優先順位で解釈する。

1. `first_fail_category_nonbody == "hook"` なら、basal hook gate が non-body 側 failure を検出したと解釈する。
2. `body_shape_pass` と body 指標は、最終 step 時点で body 変形も併発しているかを示す並列診断として sweep summary に残す。
3. `first_fail_category_nonbody == "none"` かつ body gate が fail の場合、統合 `first_fail_category` は body category とする。

## 再同定した代表条件

条件:

- `stub_mode=minimal_basal_stub`
- `body_stiffness_scale=50`
- `duration_s=0.02`
- `target=local_hook_scale`
- `values=1,8`
- 決定論的条件、Brownian off

| torque | local_hook_scale | 期待結果 | first-fail category |
| ---: | ---: | --- | --- |
| `1.2e-21 N m` | `1` | pass | `none` |
| `1.2e-21 N m` | `8` | pass | `none` |
| `4.0e-21 N m` | `1` | fail | `hook` |
| `4.0e-21 N m` | `8` | fail | `hook` |

この最小 sweep における現行境界は torque-dominated である。つまり、検証した両方の hook scale で、`1.2e-21 N m` が pass representative、`4.0e-21 N m` が hook first-fail representative である。

## モデリング上の注意

hook potential の threshold は `90 deg` のままである。現行診断 baseline の failure は hook angle error 単独ではなく、hook length / attach-first relative error が支配的である。したがって、現行モデルでは `local_hook_scale` が failure を単純に単調制御するわけではない。この sweep は failure mode を再現可能に追跡するためのものであり、最終的な物理最適値を主張するものではない。

## 検証

自動検証:

- `tests/test_motor_scale_sweep.py::test_phase24_local_hook_scale_sweep_reproduces_hook_gate`
- `tests/test_simulation.py::test_phaseb1_known_failure_case_detected_minimal_stub_005s`
- `tests/test_simulation.py::test_phaseb2_break_torque_observed_at_1_2e21_minimal_stub_005s`
- `tests/test_plot_motor_scale_collapse_heatmap.py`
