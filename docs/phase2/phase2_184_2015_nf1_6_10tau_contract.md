# Phase 2 Issue #184: 2015 project nf1–6 10τ stability screen

2015 projectの`2.5e-20 N m` per-flagellum torqueにおける、長時間の形状・遊泳診断である。
dataset採択、canonical torque選定、supported profile昇格、Phase 3 handoffは行わない。

| item | value |
| --- | --- |
| profile | `conf/sim_swim_2015.yaml` project |
| torque / scale policy | `2.5e-20 N m`; motor = reference = force torque; `reference_torque` |
| conditions | `n_flagella=1,2,3,4,5,6`; attach/phase seed `0` |
| integration | `dt_star=1e-5`, `duration_tau=10` |
| motion | RUN fixed、switchingなし、Brownian OFF |
| execution | cs10、`cs10_qualified` 3 worker、6 isolated shards |

開始前に、再解析済み#61のnf3・`2.5e-20 N m`・1τ decisionが`status=pass`でなければならない。
jobはこのdecision JSONをpreflightとして読み、PASS以外ならsimulationを起動しない。

各conditionについてstrict QC、最初のfailure criterion / 時刻 / step、wall time、steps/s、
body/flagella motion、3d+2d replayを保存する。見積りは約5.4日/condition、3 workerの二波で約11日である。
