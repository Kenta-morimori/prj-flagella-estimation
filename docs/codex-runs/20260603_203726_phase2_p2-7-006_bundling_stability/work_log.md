# Phase 2.7 bundling stability diagnostics

## 実装

- `output_sampling.fps_out_3d` を追加し、3D動画を `out_all_steps_3d=false` で間引けるようにした。
- `flagella.initial_orientation_mode=posterior_aligned` を追加し、hookは側面付着のまま、べん毛初期接線を菌体後方へ向けられるようにした。
- `step_summary.csv` に束化候補指標と菌体遊泳指標を追加した。
- `scripts/01_simulate_swimming/run_phase2_7_bundling_sweep.py` を追加し、Phase 2.7のtorque sweep summaryを出力できるようにした。

## 診断結果

短時間 `duration_s=0.02`, `n_flagella=3`, `time.dt_star=1.0e-4` では、`side_attach` / `posterior_aligned` と `1.0e-20, 2.0e-20, 3.0e-20 N m` の6条件すべてで `shape_pass_nonbody=True` だった。

0.5秒代表では、`posterior_aligned` の `0.5e-20, 1.0e-20, 2.0e-20 N m` がいずれも最終的に `hook` first-fail となった。`motor.local_spring_scale=1.2` はhook破綻を遅らせたが、`0.5e-20 N m` でも最終stepは `shape_pass_nonbody=False` だった。

後方へ向けること自体はできているが、`bundle_participation_ratio=0.0` のままであり、後方束化は未確認。

## 出力

- `outputs/phase2_7_bundling_sweep/2026-06-03_n3_dt1e4_dur0p02/phase2_7_bundling_sweep_summary.csv`
- `outputs/phase2_7_bundling_sweep/2026-06-03_n3_posterior_torque0p5e20_spring1p2_dur0p5/phase2_7_bundling_sweep_summary.csv`
- `outputs/2026-06-03/203936/render/swim3d.mp4`

## 判定

Review result は `FAIL`。実装基盤と診断は進んだが、Phase 2.7完了条件である「0.5秒以上のshape gate PASSかつ後方束化候補」は未達。
