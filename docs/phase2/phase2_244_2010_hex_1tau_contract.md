# Issue #244: 2010 hex project 1τ torque screen contract

## Scope

`2010_hex_project` は2010 projectのべん毛条件を維持し、菌体だけを六角柱30 beadsへ置換するpending evaluation candidateである。canonical model・dataset・ML policyは変更しない。

## Fixed model conditions

- profile: `year=2010`, `variant=hex_project`, `resolution=hybrid`, `implementation_status=pending`
- body: 六角柱5層、30 beads、diagonal brace ON
- flagellum: 11 beads、`ds=0.58 b`、長さ`5.8 b`、2010 projectのhelix/potential/hook/motor transmission
- attachment: `seeded_balanced_center_layer`。中心環6 slotを`attach_seed`で回転し、`n=1..6`のみ許可する。`n=4`はslot `[0,1,3,4]`でgap `[1,2,1,2]`とする。
- `phase_seed`はattachmentを変えず、初期helix phaseだけを変える。

## 1τ screen

`conf/phase2_multi_run/2010_hex_project_torque_1tau_issue244.yaml`は、`n=1,4`と motor torque `1.0, 2.0, 2.5, 3.0, 3.5 × 10^-20 N m`の10独立conditionを定義する。

- `duration_tau=1`、`dt_star=1e-4`、compact output、Brownian/switching OFF
- reference torqueは`2.5e-20 N m`に固定し、motor torqueだけを変える
- `conf/phase2_parallel/issue244_2010_hex_torque_1tau/job.yaml`は`cs10_qualified`、3 workers、全10条件のgeometry preflightを必須とする
- cs10 parallel jobは再現可能なheavy-run経路として維持する。2026-09-19のscreenはユーザー判断によりローカル直列`run_multi_run.py`で実行した。

## Evidence and next step

各condition manifestは観測されたtotal/body/flagella beads、spring segment数、segment-repulsion pair数を記録し、wall timeとsteps/sの比較に使う。

2026-09-19のlocal screen（`outputs/2026-09-19/142000/`）は10/10 conditionが1τ（10,000 steps）を完走し、各conditionに41-state archiveとtrajectoryを保存した。最終`shape_pass_nonbody`は全conditionで`True`だった一方、全conditionが最初の内部step（`4e-6 s`）で`hook` first-failを記録した。よってこのscreenはtorqueを採択せず、54-condition seed gridへ進めない診断結果とする。3D/2D grid replayは`analysis/replay/`に保存した。600τ campaignはIssue #245の範囲である。
