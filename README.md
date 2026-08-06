# prj-flagella-estimation

顕微鏡動画から細菌個体ごとのべん毛本数を推定するための開発リポジトリです。

実顕微鏡動画では個体ごとのべん毛本数ラベルを直接得にくいため、Phase 2 では 3D 物理シミュレーションと 2D 擬似顕微鏡動画生成を整備し、後続の検出・学習パイプラインへ渡せる教師データを作ります。

## 現在地

現在の主対象は **Phase 2** です。

- Phase 2: 3D 物理シミュレーションと 2D 擬似顕微鏡動画生成
- Phase 2.8: RUN 固定条件で、べん毛本数差による遊泳挙動を dataset として分析
- Phase 3 以降: 菌体検出、個体クリップ生成、べん毛本数推定モデルの学習・評価

詳細な現在地は，各Phaseのcurrent文書を参照してください。

- Phase 2: `docs/phase2/phase2_current.md`
- Phase 3: `docs/phase3/phase3_current.md`
- Phase 4: `docs/phase4/phase4_current.md`

Phase間の依存関係と作業順序は，GitHub Issues / Projectsを正本とします。

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
  config=conf/phase2_multi_run/flagella_count_behavior_v0.yaml
uv run python scripts/02_phase2_analysis/build_dataset.py \
  config=conf/phase2_multi_run/flagella_count_behavior_v0.yaml
```

CLI の詳細、override 例、後から 3D / 2D render を再生成する方法は `scripts/README.md` を参照してください。

## 主要ディレクトリ

- `conf/`: 実行設定。Phase 2 simulation のdefaultは `conf/sim_swim_2010.yaml`、複数条件 run / heatmap / replay / dataset 作成は `conf/phase2_multi_run/` を使います。

Phase 2のmodel profileは次の4つです。`2015` paper profileはmotor dynamicsを
`paper_inspired_approximation` として実装済みです。refined 120-bead geometryも実装済みで、
両2015 profileは#168の安定性評価用に実行できますが、採否未確定の`pending` profileです。

| config | role |
| --- | --- |
| `conf/sim_swim_2010.yaml` | 現行2010 project default |
| `conf/sim_swim_2010_paper.yaml` | Watari & Larson (2010) paper参照条件 |
| `conf/sim_swim_2015.yaml` | 2015 refined modelのproject採用候補 |
| `conf/sim_swim_2015_paper.yaml` | Kong et al. (2015) paper参照条件 |

2015 paper motor dynamicsは
`motor.force_distribution=hook_coupled_body_reaction` を使います。hook基部近傍の
flagellumへ駆動torqueを与え、attach beadとbody one-ringへ逆向きのzero-net-force
torque coupleを与えます。局所supportが縮退した場合だけ全body beadsへfallbackし、
support bead数とfallback使用有無をmanifestへ保存します。これは論文完全再現では
ありません。2015 project profileは比較用に既存の`root_torque_segment_couples`を
維持します。parameterとprovenanceの対応は
`docs/phase2/phase2_167_2015_paper_conditions.md`を参照してください。

2015 profilesはmanifest上で`geometry=implemented`, `simulation=evaluation_ready`を記録し、
#168向けの評価実行を許可します。`implementation_status=supported`への昇格は#168まで行いません。
project profileは3本の初期helix axisを菌体後方へ揃え、paper profileは側方初期配置を維持します。
両profileの2D確認像はPhase 3 common clipと同じ`body_capsule_orthographic_v1`を使い、
3D確認像には実時間`t`と時間scale`tau_s`を表示します。
現在の状態とmanifest provenanceはsimulationを開始せず次のように確認できます。

```bash
uv run python -c "import yaml; from pathlib import Path; from sim_swim.sim.params import SimulationConfig; cfg=SimulationConfig.from_dict(yaml.safe_load(Path('conf/sim_swim_2015_paper.yaml').read_text())); print(cfg.implementation_manifest())"
```

Phase 2 の時間指定は `time.duration.value` / `time.duration.unit` と
`time.integration.dt_star` が正本です。`duration.unit` は `s` または `tau`
を受け付けます。`time.duration_s`、`time.dt_star`、`time.dt_s` は既存互換の
deprecated key です。`time.dt_s` は出力・記録間隔であり、内部積分刻みは
`time.integration.dt_star` から決まります。

内部step数は `total_steps = ceil(duration_tau / dt_star)` です。
`step_summary.csv` は各stepの開始時刻を `total_steps` 行記録し、state列は
初期状態を含む `total_steps + 1` 件になります。そのため `ceil` により最終
state時刻が指定durationをわずかに超える場合があります。
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

- `scripts/README.md`: CLI一覧，実行例，override指定方法
- `docs/phase2/phase2_current.md`: Phase 2の現在地
- `docs/phase2/phase2_tasks.md`: Phase 2の採択判断と根拠
- `docs/codex/codex_workflow.md`: Codex作業・review result・commit / push方針
