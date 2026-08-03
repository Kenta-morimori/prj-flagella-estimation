# Phase 2 Issue #168: 2015 refined model Stage A検証契約

## 現在地

2015 project / paper profileはgeometryとdynamicsが実装済みで、Stage A評価を実行できる。
profile自体は`pending`のままとする。2026-08-03にmotor-off pilotを完了し、閾値を固定した。
現在はmotor-on RUNのユーザー実行待ちである。

Stage Aは次の順に進める。

1. brace OFFでproject/paperのmotor-off `0.1 tau`を実行する。
2. pilot結果から共通閾値案を生成し、値を固定する。
3. 固定済み閾値でbrace OFFのmotor-on RUN `1 tau`を実行する。
4. body gateが失敗したprofileだけbrace ONを比較する。
5. 自動判定後、project/paperのmotor-on replayをユーザーが目視する。

## 実行条件

| stage | duration | `dt_star` | steps | motor |
| --- | ---: | ---: | ---: | --- |
| motor-off pilot | `0.1 tau` | `1e-5` | 10,000 | forceのみOFF |
| motor-on RUN | `1 tau` | `1e-5` | 100,000 | ON |

両stageともRUN固定、polymorph switching / reversalなし、Brownian OFFとする。motor-offでも
`reference_torque_Nm=1.2e-18`を時間換算とpotential scale用に維持する。

判定用`step_summary.csv`とaggregate body diagnosticsは全step保存する。trajectoryと
`state_archive.npz`はinitial/finalを含む201 stateへ間引き、project/paper比較replayを
再シミュレーションなしで生成できるようにする。per-constraint body local CSVは全step保存しない。

2026-08-03のpaper motor-on 1-step `cProfile`では、profiler込み`0.500 s/step`のうち
segment repulsionがcumulative 70.6%、torsion有限差分が14.4%、RPY mobilityが11.0%だった。
主要ボトルネックは全segment pairを評価するspring-spring repulsionである。これは短時間の
実装診断値であり、正式なsteps/sと総時間は各user runの`performance.json`を正本とする。

## 閾値固定規則

motor-offのproject/paper最大値に対し、共通閾値案を
`min(cap, max(floor, 5 * pilot_max))`で生成する。pilot値がcapを超えた場合は閾値を緩めず、
motor-off failureとして診断する。floor/capは
`src/sim_swim/analysis/stage_a_2015_analysis.py`の`THRESHOLD_POLICY`を正本とする。

主なcap:

- body / flag / Hook spring相対誤差: `0.08`
- bend / torsion / Hook角誤差: `30 / 60 / 30 deg`
- body length / width / cross-section area変化: `10 / 10 / 20%`
- helix radius絶対誤差: `0.035 b`
- helix pitch相対誤差: `5%`
- motor body/flag action-reaction residual: `1e-8`

swimmer全体のnet force / torqueとnormalized residualは全step記録する。ただし既存のtorsion
有限差分等を含む総forceは、1-step確認でもmotor-offでforce約`3.4e-2`、torque約`4.7e-2`の
normalized residualを持つため、Stage Aで新たなゼロ残差gateは導入しない。motor-off pilotの
時系列とproject/paper差を結果へ保存し、物理的な釣り合い改善が必要なら別taskで扱う。

motor-onの各flagellumはbody-relative回転が
`max(0.01 rev, 10 * motor-off drift)`以上であることを要求する。motor body/flag torqueの
action-reaction residualは`1e-8`以下とする。131 Hz較正はStage A対象外である。

具体値はmotor-off pilot後に
`conf/phase2_validation/2015_stage_a_thresholds.yaml`へ`status: locked`として記録する。
locked前のmotor-on解析は拒否する。

## Motor-off pilot結果

`outputs/2026-08-03/111818`でproject/paperとも10,000 stepsを例外・非有限値なしで完走した。
実測はproject `3.780 steps/s`、paper `3.784 steps/s`で、各条件の所要時間は約44分だった。
両profileともbody gateはPASSし、brace ON比較は不要と判断した。

閾値案は`outputs/2026-08-03/125815/threshold_proposal.json`へ保存した。pilot値は全cap内で、
提案値を`conf/phase2_validation/2015_stage_a_thresholds.yaml`へ固定した。paper motor-offの
body-relative回転最大値は`0.8133947 rev / 0.1 tau`であり、既定のduration換算規則により
motor-on最低回転を`8.133947 rev / 1 tau`とした。この値は131 Hz較正ではなく、motor-off driftを
下回らないためのStage A検出閾値である。

## ユーザー実行

最初にprofileと条件だけを確認する。

```bash
uv run python scripts/01_simulate_swimming/run_sweep.py \
  config=conf/phase2_sweeps/2015_stage_a_motor_off.yaml \
  dry_run=true
```

motor-off pilot:

```bash
uv run python scripts/01_simulate_swimming/run_sweep.py \
  config=conf/phase2_sweeps/2015_stage_a_motor_off.yaml
```

出力されたrun rootを解析する。

```bash
uv run python scripts/01_simulate_swimming/analyze_2015_stage_a.py \
  --motor-off-run <motor-off-run-root>
```

Codexが`threshold_proposal.json` / `threshold_proposal.md`を確認し、閾値を固定した後にだけ
motor-onを実行する。

```bash
uv run python scripts/01_simulate_swimming/run_sweep.py \
  config=conf/phase2_sweeps/2015_stage_a_motor_on.yaml
```

motor-on解析:

```bash
uv run python scripts/01_simulate_swimming/analyze_2015_stage_a.py \
  --motor-off-run <motor-off-run-root> \
  --motor-on-run <motor-on-run-root> \
  --threshold-contract conf/phase2_validation/2015_stage_a_thresholds.yaml
```

replay:

```bash
uv run python scripts/01_simulate_swimming/render_shape_stability_grid_replay.py \
  input_dir=<motor-on-run-root> \
  output_dir=<motor-on-run-root>/replay \
  mode=render-only
```

## 判定とbrace分岐

body gateが失敗したprofileだけ、同じstage/profileに
`diagonal_braces_enabled=true`を指定して再実行する。nonbody failureだけではbrace比較を行わない。
brace採用、profileの`supported`昇格、Stage B移行は自動化せず、結果とユーザー目視後に判断する。

目視項目はcollapse / fly-away、body変形、Hook長・角度、helix形状、flagellumの軸中心回転、
body反作用である。目視未完了の間はIssue #168のreview resultを`FAIL`に保つ。
