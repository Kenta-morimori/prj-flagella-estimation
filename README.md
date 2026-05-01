# プロジェクト概要
顕微鏡動画からバクテリアのべん毛本数を推定するための開発リポジトリです。

## セットアップ
- Python 3.11 推奨
- 依存インストール: `uv sync`
- Git hook 有効化: `./scripts/setup_git_hooks.sh`

## pre-commit hook
- 本リポジトリは `.githooks/pre-commit` で以下を実行します。
  - `ruff format --check`
  - `ruff check`
  - `pytest`（`tests/test_*.py` が存在する場合）
- 失敗時は `git commit` をブロックします。
- 解除方法: `git config --unset core.hooksPath`

## Phase2 実行
`uv run python -m scripts.01_simulate_swimming --config conf/sim_swim.yaml`

実行後、`outputs/YYYY-MM-DD/HHMMSS/` 配下に `run.log` と `manifest.json` を含む成果物が生成されます。

## Phase2 運用ポリシー（Issue #37）
- PhaseA（`motor.torque_Nm=0`）は CSV による判定を行い、図の作成は実施しない。
- 図の作成は PhaseB1（`motor.torque_Nm!=0`）以降で開始する。
- 合否判定は `pos_all_finite` 単独ではなく、shape 系列（`shape_pass_nonbody` / body diagnostics / fail category）を優先する。

### PhaseB1 例: torque x scale 観測CSV
```bash
uv run python -m scripts.run_motor_scale_sweep \
  --target local_hook_scale \
  --values 8,2,0 \
  --torques 1e-21,4e-21,1e-20 \
  --duration 0.05 \
  --output-dir outputs/phaseb1_torque_scale
```

この実行は図を生成せず、`*_sweep_summary.csv` に `torque_Nm`, `shape_pass`, `first_fail_category` を含む観測結果を保存します。

PhaseB1 以降では、同CSVから崩壊の分布を可視化できる。

```bash
uv run python -m scripts.plot_motor_scale_collapse_heatmap \
  --summary-csv outputs/phaseb1_torque_scale/local_hook_scale_sweep_summary.csv \
  --output-dir outputs/phaseb1_torque_scale/plots
```

この可視化は `*_category_heatmap.png` に加えて、body / hook / flagella の pass/fail heatmap を出力する。
出力名はそれぞれ `*_body_pass_fail_heatmap.png`、`*_hook_pass_fail_heatmap.png`、`*_flagella_pass_fail_heatmap.png` である。
body と hook は常時出力し、flagella は flag 系の first-fail がある場合にのみ出力する。
また、正規化した `*_collapse_matrix.csv` も出力する。
