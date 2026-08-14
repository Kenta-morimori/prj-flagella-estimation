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

## Interrupted campaign の定性replay

手動停止などでconditionが未完了の場合、通常の集計は失敗する。未完了conditionを暗黙に落とさないためである。既存archiveを用いて完了conditionだけを定性確認する場合に限り、明示opt-inで次を実行する。

```bash
uv run python scripts/02_phase2_analysis/analyze_2010_torque_fixed_real_time_performance.py \
  --run-dir <campaign-root> \
  --allow-incomplete \
  --render-qualitative-replay
```

これはsimulationを再実行しない。`analysis/fixed_real_time_qualitative_replay/`に3D grid MP4と最終frame PNG、manifestを保存する。rootの`qc_summary.json`とreplay manifestには、除外したcondition ID・停止状態・不足artifactを記録する。

replayは同一実時間`0.5 s`の比較であり、各torqueの`duration_tau`は異なる。したがって、未完了conditionをPASS/FAILとして扱わず、無次元相似性、`dt_star`採択、torque採択、dataset採択の根拠には用いない。

## 2026-08-14 body stability intervention screen の解釈

`T=1e-19 N m`、`dt_star=1e-3`、`3.2 tau`で、diagonal braceの有無と
`stiffness_scales.body=[1.0, 1.5, 2.0]`を比較した。全6条件がbody shape gateを通過した。
braceを残したままbody scaleを上げると、screen内の最大spring伸びは
`0.0812 -> 0.0767 -> 0.0732`、最大body bend誤差は
`23.0 -> 20.5 -> 19.0 deg`に低下した。一方、braceを外すと同一scaleで
centerline偏差が大きく、除去を対策候補として支持する証拠は得られなかった。

ただし、このscreenのbaseline（brace ON, body scale 1.0）は、先行した
fixed-real-time runで`2.725 tau`に観測したbody spring failureを再現しなかった。
同一PCで`output.policy=debug`と`compact`を比較したところ、3.2 tauまでの
state archiveは一致したため、出力・集計形式が原因ではない。先行runとの初期形状差は
丸め誤差レベルであり、`2.5--2.8 tau`付近で増幅される数値的に敏感な領域であることが
示唆される。

従って、このscreenは「brace除去よりbody scale増加のほうが短時間指標を悪化させない」
という候補選別に限って使う。body scaleのdefault採用、torque採択、dataset採択はしない。
それらには、同一実行環境でbaselineと候補を並べた固定実時間の長時間検証と、初期位相・
attachment seedまたは実行環境を変えた再現性評価が必要である。

### 次の長時間検証（user-executed）

`conf/phase2_multi_run/2010_project_body_stability_robustness_0p5s.yaml`は、
braceを維持したbaseline（body scale `1.0`）と短時間screen候補（`2.0`）を、
phase seed `0, 1, 2`で比較する6 conditionのcampaignである。すべて
`T=1e-19 N m`, `dt_star=1e-3`, `0.5 s = 50 tau`、各50,001 stepである。
同一PCで実行し、各seedでbaseline/candidateを比較する。candidateをdefaultへ
採用するには、baselineよりfail countが増えず、body gateとmotor residualを満たし、
replay上で明瞭な非物理的変形がないことを要求する。

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/2010_project_body_stability_robustness_0p5s.yaml \
  dry_run=true
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/2010_project_body_stability_robustness_0p5s.yaml
uv run python scripts/02_phase2_analysis/render_phase2_replay.py \
  config=conf/phase2_multi_run/2010_project_body_stability_robustness_0p5s.yaml \
  run_dir=<campaign-root> \
  --mode both
```

解析はsimulationを再実行しない。`<campaign-root>`はrun出力の
`outputs/YYYY-MM-DD/HHMMSS/`であり、`analysis/body_stability_robustness_replay/`に
feature metrics CSV/PNGと3D grid replayを保存する。
