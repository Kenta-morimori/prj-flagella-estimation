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

正本は`n=1..6`、motor torque
`1.0, 2.0, 2.5, 3.0, 3.5 × 10^-20 N m`、`dt_star=1e-4,1e-3`の60 conditionとする。
既存local run `outputs/2026-09-19/142000/`と`outputs/2026-09-20/000940/`の`n=1,4`・20 conditionを再利用し、
`conf/phase2_multi_run/2010_hex_project_torque_dt_1tau_issue244.yaml`で不足40 conditionだけを追加する。

- `duration_tau=1`、compact output、Brownian/switching OFF
- reference torqueは`2.5e-20 N m`に固定し、motor torqueだけを変える
- 新規40 conditionはユーザーがlocal直列`run_multi_run.py`で実行する。既存・新規とも同一Macのため、physical QCに加えwall timeとsteps/sも比較する。
- `scripts/03_dataset_building/analyze_dataset.py --analysis-kind model-development-evaluation`は複数runを統合し、`n=1..6`別の`dt_star`（横）× torque（縦）heatmapを作る。
- hook angleはdiagnostic-onlyとし、hook length、finite、body、flag、motor QCはPASS/FAILを維持する。warning状態は持たない。

不足40 conditionは以下で起動する。

```bash
.venv/bin/python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/2010_hex_project_torque_dt_1tau_issue244.yaml \
  'sweep.include_condition_ids=[nf02__tq1p0e20__dt1e4,nf02__tq1p0e20__dt1e3,nf02__tq2p0e20__dt1e4,nf02__tq2p0e20__dt1e3,nf02__tq2p5e20__dt1e4,nf02__tq2p5e20__dt1e3,nf02__tq3p0e20__dt1e4,nf02__tq3p0e20__dt1e3,nf02__tq3p5e20__dt1e4,nf02__tq3p5e20__dt1e3,nf03__tq1p0e20__dt1e4,nf03__tq1p0e20__dt1e3,nf03__tq2p0e20__dt1e4,nf03__tq2p0e20__dt1e3,nf03__tq2p5e20__dt1e4,nf03__tq2p5e20__dt1e3,nf03__tq3p0e20__dt1e4,nf03__tq3p0e20__dt1e3,nf03__tq3p5e20__dt1e4,nf03__tq3p5e20__dt1e3,nf05__tq1p0e20__dt1e4,nf05__tq1p0e20__dt1e3,nf05__tq2p0e20__dt1e4,nf05__tq2p0e20__dt1e3,nf05__tq2p5e20__dt1e4,nf05__tq2p5e20__dt1e3,nf05__tq3p0e20__dt1e4,nf05__tq3p0e20__dt1e3,nf05__tq3p5e20__dt1e4,nf05__tq3p5e20__dt1e3,nf06__tq1p0e20__dt1e4,nf06__tq1p0e20__dt1e3,nf06__tq2p0e20__dt1e4,nf06__tq2p0e20__dt1e3,nf06__tq2p5e20__dt1e4,nf06__tq2p5e20__dt1e3,nf06__tq3p0e20__dt1e4,nf06__tq3p0e20__dt1e3,nf06__tq3p5e20__dt1e4,nf06__tq3p5e20__dt1e3]'
```

## Stage 2: seed grid and Δt convergence

Stage 1の60 cell heatmapが全てPASSであることを確認した後に限り、
`conf/phase2_multi_run/2010_hex_project_seed_grid_1tau_issue244.yaml`で
`T=2.5e-20 N m/flagellum`・`dt_star=1e-4`、`n=1..6`、attach / phase seed各0..2の
54 conditionを実行する。

54 conditionが通過した後に限り、
`conf/phase2_multi_run/2010_hex_project_dt_convergence_1tau_issue244.yaml`で
`n=3,6`・seed 0・`dt_star=1e-4,5e-5`の4 conditionを比較する。Issue #244の全計画は
Stage 1の60 + seed grid 54 + convergence 4 = 118 conditionとする。

## Evidence and next step

各condition manifestは観測されたtotal/body/flagella beads、spring segment数、segment-repulsion pair数を記録し、wall timeとsteps/sの比較に使う。

2026-09-19のlocal baseline（`outputs/2026-09-19/142000/`）は10/10 conditionが1τ（10,000 steps）を完走し、各conditionに41-state archiveとtrajectoryを保存した。最終`shape_pass_nonbody`は全conditionで`True`だった一方、全conditionが最初の内部step（`4e-6 s`）で`hook` first-failを記録した。hook angleは後方束化や長さ破綻を意味しないため、診断として保持しつつPASS/FAILから除外する。600τ campaignはIssue #245の範囲である。

2026-09-20のlocal `dt_star=1e-3` 10 condition（`outputs/2026-09-20/000940/`）は、全conditionが1τ（1,000 steps）を完走し、41-state archive、trajectory、3D/2D replayを保存した。共通評価器は既存20 conditionと新規40 conditionを`model_development_evaluation/`へ統合し、PNG heatmapと固定camera・各べん毛軸付きのΔt別MP4を出力する。

Stage 2はcs10 user-run対象であり、ユーザーが明示的に開始を許可した後に限る。54 condition jobは`conf/phase2_parallel/issue244_2010_hex_seed_grid_1tau/job.yaml`、その後の4 condition刻み比較jobは`conf/phase2_parallel/issue244_2010_hex_dt_convergence_1tau/job.yaml`を使用する。いずれも3 workers、全condition geometry preflightを必須とする。
