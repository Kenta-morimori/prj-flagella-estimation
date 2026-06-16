# Phase 2.7 axis-alignment bundling definition

## 実施内容

- Phase 2.7 の束化定義を、近接束ではなく複数べん毛の螺旋中心軸方向の安定整列へ更新した。
- `step_summary.csv` に軸整列指標を追加した。
- `run_phase2_7_bundling_sweep.py` の分類を `axis_aligned_stable`, `hook_wrapped_axis_aligned`, `axis_not_aligned`, `collapse` に更新した。
- 各run directoryに `flag_helix_axis_angles_timeseries.png` を自動生成するようにした。
- docs/phase2 のタスク定義と診断docを更新した。

## 代表結果

条件:

- `n_flagella=3`
- `motor.torque_Nm=2.5e-20`
- `duration_s=0.5`
- `time.dt_star=1.0e-4`
- `flagella.initial_helix_axis_from_rear_deg=0`

結果:

- `phase27_class=hook_wrapped_axis_aligned`
- `phase27_axis_alignment_stable=True`
- `phase27_axis_alignment_stable_fraction=1.0`
- `flag_helix_axis_mean_deviation_deg_max=11.8841`
- `first_fail_category_nonbody=hook`

出力:

- `outputs/phase2_7_axis_alignment_v1/phase2_7_bundling_sweep_summary.csv`
- `outputs/phase2_7_axis_alignment_v1/helix_axis_angle_0deg/n_3/torque_2.50e-20/flag_helix_axis_angles_timeseries.png`

## 検証

- `uv run ruff check ...`
- `uv run ruff format --check ...`
- `uv run pytest tests/test_helix_axis.py tests/test_phase2_7_bundling_sweep.py tests/test_simulation.py -q`
- short sweep smoke
- representative 0.5s sweep
- `git diff --check`
