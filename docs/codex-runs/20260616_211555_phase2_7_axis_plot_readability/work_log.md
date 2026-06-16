# Phase 2.7 axis plot readability

## 実施内容

- `flag_helix_axis_angles_timeseries.png` の上段を `0-90 deg` 固定にした。
- 上段の `rear axis (0 deg)` と `side axis (90 deg)` にlegendを追加した。
- 下段の赤点線に `alignment threshold (15 deg)` のlegendを追加した。
- 代表条件のplotを再生成した。

## 検証

- `uv run ruff check scripts/01_simulate_swimming/run_phase2_7_bundling_sweep.py tests/test_phase2_7_bundling_sweep.py`
- `uv run ruff format --check scripts/01_simulate_swimming/run_phase2_7_bundling_sweep.py tests/test_phase2_7_bundling_sweep.py`
- `uv run pytest tests/test_phase2_7_bundling_sweep.py -q`
- 代表条件sweepの再実行
