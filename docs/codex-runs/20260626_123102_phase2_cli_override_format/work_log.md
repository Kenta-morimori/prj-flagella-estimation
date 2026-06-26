# Phase 2 CLI override format

## 目的

`scripts.01_simulate_swimming` の user-facing command examples と manifest 記録を `KEY=VALUE` override 形式へ寄せ，`--duration-s` / `--fps-out` と `time.duration_s=...` の混在を避ける。

## 実施内容

- `--duration-s` / `--fps-out` の help text を legacy shorthand と明記し，`time.duration_s=...` / `output_sampling.fps_out_2d=...` を推奨する文面へ変更した。
- `manifest.json` の `input` に，生の `cli_overrides`，legacy shorthand 展開後の `legacy_shorthand_overrides`，実効値の `effective_overrides` を保存するようにした。
- Typer decorated 関数をテストから直接呼ぶ場合の `OptionInfo` default を未指定として正規化した。
- `scripts/README.md`，`AGENTS.md`，`docs/codex/codex_workflow.md`，`docs/phase2/phase2_current.md`，`docs/phase2/phase2_6_helix_retention_gate.md` の単発 simulation 例・運用メモを override 形式へ更新した。

## 検証

- `uv run ruff format --check scripts/01_simulate_swimming/01_simulate_swimming.py tests/test_e2e_cli.py`
- `uv run ruff check scripts/01_simulate_swimming/01_simulate_swimming.py tests/test_e2e_cli.py`
- `uv run pytest tests/test_e2e_cli.py -q`
- commit hook:
  - `uv run ruff format --check .`
  - `uv run ruff check .`
  - `uv run pytest -q`
- `uv run python scripts/01_simulate_swimming/01_simulate_swimming.py --help`
- `uv run python scripts/01_simulate_swimming/01_simulate_swimming.py -- time.duration_s=0.001 time.dt_star=1.0e-4 output.base_dir=/private/tmp/phase2_duration_override_smoke`

Smoke output:

- `/private/tmp/phase2_duration_override_smoke/2026-06-26/123037/manifest.json`
- `/private/tmp/phase2_duration_override_smoke/2026-06-26/123037/run.log`

確認結果:

- manifest の `input.overrides` に `time.duration_s=0.001`，`time.dt_star=1.0e-4`，`output.base_dir=/private/tmp/phase2_duration_override_smoke` が保存された。
- manifest の `input.effective_overrides.time.duration_s` は `0.001`，`input.effective_overrides.time.dt_star` は `0.0001`。
- run.log に `Overrides: {'time': {'duration_s': 0.001, 'dt_star': 0.0001}, ...}`，`dt_star=1.000000e-04`，`duration_s=1.000000e-03` が出力された。

## レビュー判断

CLI・manifest・文書例の整理であり，物理モデルや描画内容は変更していない。ユーザー目視レビューと ADR は不要。
