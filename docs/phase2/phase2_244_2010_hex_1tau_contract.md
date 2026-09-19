# Issue #244: 2010 hex project 1τ torque screen contract

## Scope

`2010_hex_project` は2010 projectのべん毛条件を維持し、菌体だけを六角柱30 beadsへ置換するpending evaluation candidateである。canonical model・dataset・ML policyは変更しない。

## Fixed model conditions

- profile: `year=2010`, `variant=hex_project`, `resolution=hybrid`, `implementation_status=pending`
- body: 六角柱5層、30 beads、diagonal brace ON
- flagellum: 11 beads、`ds=0.58 b`、長さ`5.8 b`、2010 projectのhelix/potential/hook/motor transmission
- attachment: `seeded_balanced_center_layer`。中心環6 slotを`attach_seed`で回転し、`n=1..6`のみ許可する。`n=4`はslot `[0,1,3,4]`でgap `[1,2,1,2]`とする。
- `phase_seed`はattachmentを変えず、初期helix phaseだけを変える。

## Stage 1: 1τ torque--Δt screen

既存local run `outputs/2026-09-19/142000/` は、`n=1,4`とmotor torque
`1.0, 2.0, 2.5, 3.0, 3.5 × 10^-20 N m`、`dt_star=1e-4`の10 conditionである。
これを再実行せず、`conf/phase2_multi_run/2010_hex_project_torque_dt_1tau_issue244.yaml`
で同じ10 torque--count cellの`dt_star=1e-3`を追加する。Stage 1全体は20独立conditionである。

- `duration_tau=1`、compact output、Brownian/switching OFF
- reference torqueは`2.5e-20 N m`に固定し、motor torqueだけを変える
- 新規10 conditionはユーザーがlocal直列`run_multi_run.py`で実行する。既存・新規とも同一Macのため、physical QCに加えwall timeとsteps/sも比較する。
- `scripts/03_dataset_building/analyze_dataset.py --analysis-kind issue244-torque-dt`は両runを統合し、`n=1`・`n=4`別のtorque × `dt_star`（2×5）heatmapを作る。
- strict simulator gateは変更しない。finite/body/flag/hook length/motor diagnosticsが正常で、最初の内部stepだけのhook angle違反、最終nonbody pass、failure sample数1をすべて満たすものだけをIssue #244の解析でwarningとする。

新規 `dt_star=1e-3` の10 conditionは以下で起動する。

```bash
.venv/bin/python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/2010_hex_project_torque_dt_1tau_issue244.yaml \
  'sweep.include_condition_ids=[nf01__tq1p0e20__dt1e3,nf01__tq2p0e20__dt1e3,nf01__tq2p5e20__dt1e3,nf01__tq3p0e20__dt1e3,nf01__tq3p5e20__dt1e3,nf04__tq1p0e20__dt1e3,nf04__tq2p0e20__dt1e3,nf04__tq2p5e20__dt1e3,nf04__tq3p0e20__dt1e3,nf04__tq3p5e20__dt1e3]'
```

## Stage 2: seed grid and Δt convergence

Stage 1のheatmapでwarning以外のfailがないことを確認した後に限り、
`conf/phase2_multi_run/2010_hex_project_seed_grid_1tau_issue244.yaml`で
`T=2.5e-20 N m/flagellum`・`dt_star=1e-4`、`n=1..6`、attach / phase seed各0..2の
54 conditionを実行する。

54 conditionが通過した後に限り、
`conf/phase2_multi_run/2010_hex_project_dt_convergence_1tau_issue244.yaml`で
`n=3,6`・seed 0・`dt_star=1e-4,5e-5`の4 conditionを比較する。Issue #244の全計画は
Stage 1の20 + seed grid 54 + convergence 4 = 78 conditionとする。

## Evidence and next step

各condition manifestは観測されたtotal/body/flagella beads、spring segment数、segment-repulsion pair数を記録し、wall timeとsteps/sの比較に使う。

2026-09-19のlocal baseline（`outputs/2026-09-19/142000/`）は10/10 conditionが1τ（10,000 steps）を完走し、各conditionに41-state archiveとtrajectoryを保存した。最終`shape_pass_nonbody`は全conditionで`True`だった一方、全conditionが最初の内部step（`4e-6 s`）で`hook` first-failを記録した。Stage 1統合解析でその単発hook angle transientを明示的にwarning/failへ分類し、torque・Δtの採択根拠を残す。3D/2D grid replayは`analysis/replay/`に保存した。600τ campaignはIssue #245の範囲である。
