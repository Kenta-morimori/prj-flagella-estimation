# Phase 2.7 P2-7-006 角度sweep・反発診断ログ

## 目的

複数べん毛で、形状崩壊せず後方束化する条件が見つかるかを確認する。特に、初期接線を菌体後方軸から少し開いた条件にし、束化しない原因がべん毛同士の反発・接触なのか、束化駆動不足なのかを切り分ける。

## 実装

- `flagella.initial_tangent_vs_rear_deg` を追加した。
  - `0 deg`: 菌体後方軸。
  - `10 deg`: 主候補。後方を向きつつ少し外側へ開く。
  - `90 deg`: 従来 `side_attach` 相当。
- `run_phase2_7_bundling_sweep.py` に `--tangent-angles-deg` を追加した。
- `step_summary.csv` に flagella-flagella segment distance / close-contact count / repulsion force 診断列を追加した。
- `docs/phase2/phase2_7_bundling_stability_plan.md`, `docs/phase2/phase2_tasks.md`, `docs/PROJECT_PLAN.md` に結果を記録した。

## 実行

短時間screening:

```bash
uv run python scripts/01_simulate_swimming/run_phase2_7_bundling_sweep.py --duration-s 0.02 --dt-star 1.0e-4 --torques 0.5e-20,1.0e-20 --n-flagella 3 --tangent-angles-deg 0,10,30,60,90 --output-dir outputs/phase2_7_bundling_sweep/2026-06-04_n3_angle_screen_dur0p02 render.save_frames_3d=false render.save_frames_2d=false output_sampling.out_all_steps_3d=false output_sampling.fps_out_3d=25 motor.local_spring_scale=1.2
```

0.5秒代表条件:

```bash
uv run python scripts/01_simulate_swimming/run_phase2_7_bundling_sweep.py --duration-s 0.5 --dt-star 1.0e-4 --torques 0.5e-20 --n-flagella 3 --tangent-angles-deg 0,10,30 --output-dir outputs/phase2_7_bundling_sweep/2026-06-04_n3_angle_representative_dur0p5 render.save_frames_3d=false render.save_frames_2d=false output_sampling.out_all_steps_3d=false output_sampling.fps_out_3d=25 motor.local_spring_scale=1.2
```

10度でのトルク上昇:

```bash
uv run python scripts/01_simulate_swimming/run_phase2_7_bundling_sweep.py --duration-s 0.5 --dt-star 1.0e-4 --torques 1.0e-20,2.0e-20 --n-flagella 3 --tangent-angles-deg 10 --output-dir outputs/phase2_7_bundling_sweep/2026-06-04_n3_angle10_torque_escalation_dur0p5 render.save_frames_3d=false render.save_frames_2d=false output_sampling.out_all_steps_3d=false output_sampling.fps_out_3d=25 motor.local_spring_scale=1.2
```

10度でのhook補強:

```bash
uv run python scripts/01_simulate_swimming/run_phase2_7_bundling_sweep.py --duration-s 0.5 --dt-star 1.0e-4 --torques 1.0e-20,2.0e-20 --n-flagella 3 --tangent-angles-deg 10 --output-dir outputs/phase2_7_bundling_sweep/2026-06-04_n3_angle10_spring2_dur0p5 render.save_frames_3d=false render.save_frames_2d=false output_sampling.out_all_steps_3d=false output_sampling.fps_out_3d=25 motor.local_spring_scale=2.0
```

## 結果

- 短時間screeningでは全10条件が shape gate PASS。ただし全条件 `no_bundle`, `bundle_participation_ratio=0.0`, `flag_flag_close_pair_count=0`。
- `initial_tangent_vs_rear_deg=10`, `torque_Nm=0.5e-20`, `local_spring_scale=1.2`, `duration_s=0.5` は shape gate PASS。ただし束化なし。
- `initial_tangent_vs_rear_deg=10`, `torque_Nm=1.0e-20`, `local_spring_scale=2.0`, `duration_s=0.5` も shape gate PASS。ただし束化なし。
- `1.0e-20` を `local_spring_scale=1.2` で回す条件、および `2.0e-20` 条件では、束化より先に hook length gate が破綻した。

## 判定

今回の探索範囲では、後方束化成功条件は見つからなかった。

原因は、形状が保てるトルク帯では `no_bundle_drive`、トルクを上げた条件では `hook_drift` が先行すること。flagella-flagella repulsion は今回の失敗主因ではない。

