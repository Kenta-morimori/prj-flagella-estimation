# Phase 2 Issue #168: 2015 refined model Stage A検証契約

## 現在地

2015 project / paper profileはgeometryとdynamicsが実装済みで、Stage A評価を実行できる。
profile自体は`pending`のままとする。2026-08-03にmotor-off pilot、dt reference群、canonical
motor-on RUNを完了した。現在は、project profileにおける`dt_star=1e-4`の正式採用条件と、
torque sweepを用いた回転・遊泳安定性評価を確定する段階である。

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

各conditionには軽量な`run_summary.json`を生成する。通常の確認は`manifest.json`と
`run_summary.json`を先に読み、詳細が必要な場合だけbounded `inspect_step_summary.py`を使用する。
`step_summary.csv`とaggregate body diagnosticsは全step保存し、trajectoryと`state_archive.npz`は
initial/finalを含む201 stateへ間引いて、project/paper比較replayを再シミュレーションなしで
生成できるようにする。per-constraint body local CSVは全step保存しない。

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

motor body/flag torqueのaction-reaction residualは`1e-8`以下とする。131 Hz較正はStage A対象外である。
motor-off driftから導いた最低回転量は、motor-onの物理的な採否基準としては使用しない。motor-onの
回転gateは、net revolution、方向一貫性、shape gate、遊泳運動を分離した評価契約として、torque sweep
実行前に固定する。

具体値はmotor-off pilot後に
`conf/phase2_validation/2015_stage_a_thresholds.yaml`へ`status: locked`として記録する。
locked前のmotor-on解析は拒否する。

## Motor-off pilot結果

`outputs/2026-08-03/111818`でproject/paperとも10,000 stepsを例外・非有限値なしで完走した。
実測はproject `3.780 steps/s`、paper `3.784 steps/s`で、各条件の所要時間は約44分だった。
両profileともbody gateはPASSし、brace ON比較は不要と判断した。

閾値案は`outputs/2026-08-03/125815/threshold_proposal.json`へ保存した。pilot値は全cap内で、
提案値を`conf/phase2_validation/2015_stage_a_thresholds.yaml`へ固定した。なお、paper motor-offの
body-relative回転最大値から導いた旧`8.133947 rev / 1 tau`は、motor-onの物理採否gateとしては
撤回し、reference診断値としてのみ保持する。

## 実行済み結果と確定方針

- `dt_star=1e-4`と`1e-5`の先頭`0.1 tau`比較では、project profileは全flagellumで回転方向一致・
  回転量差10%以内、paper profileは回転量差22--34%で不合格だった。
- body姿勢差と並進除去後bead差は両profileで基準内だった。従って`dt_star=1e-4`はproject用の
  非canonical referenceとして維持するが、project/paper共通の2015 defaultには採用しない。
- canonical motor-on `1 tau`ではbody shapeは両profileで維持された。一方、pitch誤差、projectの
  torsion・motor torque balance、ならびに旧回転gateが不合格だった。旧回転gateは上記のとおり採否に
  用いない。pitch/torsion/torque balanceの挙動はtorque sweepで再評価する。
- `grid_swim3d.mp4`の2015表示は`τ`基準時間、固定小数の秒、積分step数を併記する。

## ユーザー実行

### torque screen: 広域・短時間

基準`1.2e-18 N m`の`0.1, 0.2, 0.5, 1, 2, 5, 10`倍を、project/paper両profileで
`dt_star=1e-5`、`0.02 tau`（2,000 steps）実行する。これは1--2桁の広域screenであり、
長時間の採否をここだけで決めない。有限性、body/nonbody shape、motor action-reaction、
回転・遊泳診断をconditionごとの`manifest.json`と`run_summary.json`から確認する。

```bash
uv run python scripts/01_simulate_swimming/run_sweep.py \
  config=conf/phase2_sweeps/2015_stage_a_torque_screen.yaml \
  dry_run=true

uv run python scripts/01_simulate_swimming/run_sweep.py \
  config=conf/phase2_sweeps/2015_stage_a_torque_screen.yaml
```

screen後、明白な破綻を除いた連続torque帯を`1 tau`で評価する。`dt_star=1e-4`の比較は、
この長時間評価を満たしたprofile・torqueだけで`1e-5`と対にして行う。採用はprofile・torque別に
判断し、paper profileにも同等性が確認されるまでは2015共通canonicalを`1e-5`から変更しない。

### `dt_star=1e-4`参考・採用候補評価

canonical motor-onの前に、10倍粗い内部刻みを独立した非canonical参考条件として評価する。
既存のlocked閾値は変更せず、次の順に実行する。

```bash
uv run python scripts/01_simulate_swimming/run_sweep.py \
  config=conf/phase2_sweeps/2015_stage_a_dt1e4_motor_off_reference.yaml

uv run python scripts/01_simulate_swimming/run_sweep.py \
  config=conf/phase2_sweeps/2015_stage_a_dt1e4_motor_on_reference.yaml

uv run python scripts/01_simulate_swimming/run_sweep.py \
  config=conf/phase2_sweeps/2015_stage_a_dt1e5_motor_on_short_reference.yaml
```

各run rootを指定して解析する。

```bash
uv run python scripts/01_simulate_swimming/analyze_2015_stage_a.py \
  --motor-off-run outputs/2026-08-03/111818 \
  --threshold-contract conf/phase2_validation/2015_stage_a_thresholds.yaml \
  --coarse-motor-off-run <dt1e-4-motor-off-run-root> \
  --coarse-motor-on-run <dt1e-4-motor-on-run-root> \
  --fine-motor-on-short-run <dt1e-5-motor-on-0.1tau-run-root>
```

`dt_sensitivity_decision.json` / `.md` / comparison CSVへ安定性、刻み感度、高速化率を保存する。
`reference_stable`には既存locked gateの通過を要求する。`adoption_candidate`にはさらにmotor-on
先頭`0.1 tau`で各flagellumの回転方向一致・回転量差10%以内、body姿勢差5度以内、body重心並進を
除いた全beadのRMS差`0.1 b`以内・最大差`0.25 b`以内を要求する。採用候補でもreplay目視とADRなしに
2015 defaultを変更しない。3 run群のproject/paper合計実行時間は約3時間を目安とする。

### Canonical Stage A

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
