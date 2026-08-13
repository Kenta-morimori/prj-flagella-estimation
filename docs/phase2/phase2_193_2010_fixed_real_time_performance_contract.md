# Phase 2 Issue #193: 2010 fixed-real-time performance characterization

## Purpose

このcontractは、2010 projectのexperiment-only torque-linked time-scale条件で、固定実時間`0.5 s`のstep数とsimulation-loop性能を記録する。無次元時間が同じであること、束化、または`dt_star`の採択を検証するものではない。

## Conditions

- per-flagellum torque: `1e-21, 2.5e-20, 1e-19, 1.2e-18 N m`
- `dt_star=1e-3`、Brownian OFF、`motor.body_reaction_full_vector=true`
- `duration_s=0.5`。`duration_tau`とは排他的である。
- 期待する数理上のstep数は`500, 12,500, 50,000, 600,000`である。内部の`ceil`規則により、浮動小数点境界では実manifestの`total_steps`が1 step多くなる場合がある。実測値はmanifestを正本とする。

`600 tau`を実行する場合はtorqueを小さくしてもstep数は減らない。`dt_star=1e-3`では常に600,000 stepである。小torqueのstep削減は、固定実時間`0.5 s`で到達する`duration_tau`が短くなることによるものであり、同一無次元時間の再現を意味しない。

## Execution and outputs

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/2010_project_torque_fixed_real_time_0p5s_performance.yaml \
  dry_run=true
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/2010_project_torque_fixed_real_time_0p5s_performance.yaml
uv run python scripts/02_phase2_analysis/analyze_2010_torque_fixed_real_time_performance.py \
  --run-dir <campaign-root>
```

run rootは`outputs/YYYY-MM-DD/HHMMSS/phase2_issue193_2010_torque_fixed_real_time_0p5s_performance/`である。`performance_summary.csv`はconditionごとの実時間、`duration_tau`、`tau_s`、内部刻み、step数、wall time、steps/s、必須safety statusを保存する。

全conditionで完走、finite、body/non-body shape、motor action-reaction residualを確認する。評価はperformance-onlyであり、`dt_star=1e-4` / `1e-5`との一致度、2015 project、dataset採択には使わない。
