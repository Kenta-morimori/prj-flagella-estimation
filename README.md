# prj-flagella-estimation

顕微鏡動画から細菌個体ごとのべん毛本数を推定するための開発リポジトリです。

実顕微鏡動画では個体ごとのべん毛本数ラベルを直接得にくいため、Phase 2 では 3D 物理シミュレーションと 2D 擬似顕微鏡動画生成を整備し、後続の検出・学習パイプラインへ渡せる教師データを作ります。

## 現在地

現在の主対象は **Phase 2** です。

- Phase 2: 3D 物理シミュレーションと 2D 擬似顕微鏡動画生成
- Phase 2.8: RUN 固定条件で、べん毛本数差による遊泳挙動を dataset として分析
- Phase 3 以降: 菌体検出、個体クリップ生成、べん毛本数推定モデルの学習・評価

詳細な現在地は `docs/phase2/phase2_current.md`、全体計画は `docs/PROJECT_PLAN.md` を参照してください。

## セットアップ

- Python 3.11 推奨
- 依存インストール: `uv sync`
- Git hook 有効化: `./scripts/setup_git_hooks.sh`

Git hook は commit 時に `ruff format --check .`、`ruff check .`、`pytest -q -m light` を実行します。full pytest は PR 前、CI、モデル・出力仕様・pipeline 変更時に明示実行します。hook で full pytest まで回したい場合は `FULL_TEST=1 git commit ...` を使います。解除が必要な場合は `git config --unset core.hooksPath` を使います。

## Quick Start

単発の Phase 2 simulation:

```bash
uv run python -m scripts.01_simulate_swimming
```

Phase 2.8 の診断用 dataset 作成:

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/flagella_count_behavior_diagnostic.yaml
uv run python scripts/02_phase2_analysis/build_dataset.py \
  config=conf/phase2_multi_run/flagella_count_behavior_diagnostic.yaml
```

CLI の詳細、override 例、後から 3D / 2D render を再生成する方法は `scripts/README.md` を参照してください。

## 主要ディレクトリ

- `conf/`: 実行設定。Phase 2 simulation は `conf/sim_swim.yaml`、複数条件 run / heatmap / replay / dataset 作成は `conf/phase2_multi_run/` を使います。
- `scripts/`: ユーザー向け CLI entrypoints。詳細は `scripts/README.md` を参照してください。
- `src/`: 再利用可能な実装本体。Phase 2 simulation は `src/sim_swim/` が中心です。
- `docs/phase2/`: Phase 2 の現在地、task status、設計・検証記録。
- `outputs/`: 実行成果物。基本的に Git 管理外です。
- `tests/`: pytest tests。

## 出力と再現性

通常の simulation run は `outputs/YYYY-MM-DD/HHMMSS/` 配下に `run.log`、`manifest.json`、`step_summary.csv`、render 出力などを保存します。

Phase 2.8 diagnostic dataset は標準で以下に出力します。

- `outputs/phase2_multi_run/<profile>/`
- `outputs/phase2_analysis/flagella_count_behavior/datasets/<dataset_id>/`

Phase 2.8 dataset では、実行時の設定を `campaign_config_used.yaml` と manifest の `effective_campaign_config` に残します。raw condition には `trajectory.csv` と `state_archive.npz` を保存し、後から replay render を再生成できます。

## テストと品質確認

```bash
uv run ruff format --check .
uv run ruff check .
uv run pytest -q -m light
uv run pytest -q
```

軽量 pre-commit 対象の pytest だけ確認する場合:

```bash
uv run pytest -q -m light
```

Phase 2 analysis 周辺だけを確認する場合:

```bash
uv run pytest tests/test_flagella_count_behavior_dataset.py
```

## 詳細ドキュメント

- `scripts/README.md`: CLI 一覧、実行例、override 指定方法
- `docs/phase2/phase2_current.md`: Phase 2 の現在地
- `docs/phase2/phase2_tasks.md`: Phase 2 の accepted task status
- `docs/PROJECT_PLAN.md`: プロジェクト全体の流れとフェーズ定義
- `docs/codex/codex_workflow.md`: Codex 作業・review result・commit/push 方針
