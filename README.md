# prj-flagella-estimation

顕微鏡動画から細菌個体ごとのべん毛本数を推定するための開発リポジトリです．

実顕微鏡動画では個体ごとのべん毛本数ラベルを直接得にくいため，3D物理simulationと2D擬似顕微鏡動画から教師データを生成し，後続の検出・学習pipelineへ渡します．

## 現在地

各Phaseの現在地は，次のcurrent文書を正本とします．

* Phase 2: `docs/phase2/phase2_current.md`
* Phase 3: `docs/phase3/phase3_current.md`
* Phase 4: `docs/phase4/phase4_current.md`

各Phaseで採択した判断とその根拠は，対応する`phaseX_tasks.md`へ記録します．

Phase間の依存関係，作業順序，詳細なacceptance criteriaは，GitHub Issues / Projectsを正本とします．

## セットアップ

* Python 3.11推奨
* 依存関係のinstall: `uv sync`
* Git hookの有効化: `./scripts/setup_git_hooks.sh`

Git hookはcommit時に，次の確認を実行します．

* `ruff format --check .`
* `ruff check .`
* `pytest -q -m light`

full pytestはPR前，CI，model・出力仕様・pipeline変更時に明示的に実行します．

hookでfull pytestまで実行する場合は，次のようにcommitします．

```bash
FULL_TEST=1 git commit ...
```

Git hookを解除する場合は，次を実行します．

```bash
git config --unset core.hooksPath
```

## Quick Start

単発のPhase 2 simulation:

```bash
uv run python -m scripts.01_simulate_swimming
```

Phase 2.8の診断用dataset作成:

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/flagella_count_behavior_v0.yaml

uv run python scripts/03_dataset_building/build_dataset.py \
  config=conf/phase2_multi_run/flagella_count_behavior_v0.yaml
```

CLIの詳細，override例，後から3D／2D renderを再生成する方法は，`scripts/README.md`を参照してください．

## 主要ディレクトリ

* `conf/`: 実行設定
* `scripts/`: ユーザー向けCLI entrypoints
* `src/`: 再利用可能な実装本体
* `docs/phase2/`: Phase 2の現在地，採択判断，現役contract・検証記録
* `docs/phase3/`: Phase 3の現在地，採択判断，現役contract・検証記録
* `docs/phase4/`: Phase 4の現在地，採択判断，現役contract・検証記録
* `outputs/`: 実行成果物．基本的にGit管理外
* `tests/`: pytest tests

Phase 2 simulationのdefault configは，`conf/sim_swim_2010.yaml`です．

複数条件run，heatmap，replay，dataset作成には，主に`conf/phase2_multi_run/`を使用します．

## Phase 2 model profiles

Phase 2のmodel profileは次の4つです．

`2015` paper profileでは，motor dynamicsを`paper_inspired_approximation`として実装しています．refined 120-bead geometryも実装済みです．

両2015 profileはIssue #168の安定性評価用に実行できますが，採否未確定の`pending` profileです．

| config                          | role                            |
| ------------------------------- | ------------------------------- |
| `conf/sim_swim_2010.yaml`       | 現行2010 project default          |
| `conf/sim_swim_2010_paper.yaml` | Watari & Larson（2010）のpaper参照条件 |
| `conf/sim_swim_2015.yaml`       | 2015 refined modelのproject採用候補  |
| `conf/sim_swim_2015_paper.yaml` | Kong et al.（2015）のpaper参照条件     |

2015 paper motor dynamicsでは，次の設定を使用します．

```text
motor.force_distribution=hook_coupled_body_reaction
```

hook基部近傍のflagellumへ駆動torqueを与え，attach beadとbody one-ringへ逆向きのzero-net-force torque coupleを与えます．

局所supportが縮退した場合のみ，全body beadsへfallbackします．support bead数とfallbackの使用有無はmanifestへ保存します．

この実装は論文の完全再現ではありません．

2015 project profileは，比較用として既存の`root_torque_segment_couples`を維持します．parameterとprovenanceの対応は，次を参照してください．

```text
docs/phase2/phase2_167_2015_paper_conditions.md
```

2015 profilesはmanifest上で，次の状態を記録します．

```text
geometry=implemented
simulation=evaluation_ready
```

`implementation_status=supported`への昇格は，Issue #168の評価完了後に判断します．

project profileは，3本の初期helix axisを菌体後方へ揃えます．paper profileは，側方初期配置を維持します．

両profileの2D確認像では，Phase 3 common clipと同じ`body_capsule_orthographic_v1`を使用します．3D確認像には，実時間`t`と時間scale`tau_s`を表示します．

現在の状態とmanifest provenanceは，simulationを開始せずに次のコマンドで確認できます．

```bash
uv run python -c "import yaml; from pathlib import Path; from sim_swim.sim.params import SimulationConfig; cfg=SimulationConfig.from_dict(yaml.safe_load(Path('conf/sim_swim_2015_paper.yaml').read_text())); print(cfg.implementation_manifest())"
```

## Phase 2の時間指定

Phase 2の時間指定では，次を正本とします．

```text
time.duration.value
time.duration.unit
time.integration.dt_star
```

`time.duration.unit`には，`s`または`tau`を指定できます．

次のkeyは既存configとの互換性を維持するためのdeprecated keyです．

```text
time.duration_s
time.dt_star
time.dt_s
```

`time.dt_s`は出力・記録間隔です．内部積分刻みは，`time.integration.dt_star`から決定します．

内部step数は，次の式で決定します．

```text
total_steps = ceil(duration_tau / dt_star)
```

`step_summary.csv`は，各stepの開始時刻を`total_steps`行記録します．state列には初期状態を含む`total_steps + 1`件が含まれます．

`ceil`を使用するため，最終state時刻が指定durationをわずかに超える場合があります．

## 出力と再現性

通常のsimulation runは，次のディレクトリ以下へ出力します．

```text
outputs/YYYY-MM-DD/HHMMSS/
```

主な出力には，次が含まれます．

* `run.log`
* `manifest.json`
* `step_summary.csv`
* render出力

Phase 2.8 diagnostic datasetは，標準で次へ出力します．

```text
outputs/phase2_multi_run/<profile>/
outputs/phase2_analysis/flagella_count_behavior/datasets/<dataset_id>/
```

Phase 2.8 datasetでは，実行時の設定を次へ保存します．

* `campaign_config_used.yaml`
* manifest内の`effective_campaign_config`

raw conditionには，`trajectory.csv`と`state_archive.npz`を保存します．これらを使用して，後からreplay renderを再生成できます．

## テストと品質確認

```bash
uv run ruff format --check .
uv run ruff check .
uv run pytest -q -m light
uv run pytest -q
```

軽量pre-commit対象のpytestのみ確認する場合:

```bash
uv run pytest -q -m light
```

Phase 2 analysis周辺のみ確認する場合:

```bash
uv run pytest tests/test_flagella_count_behavior_dataset.py
```

## 詳細ドキュメント

* `scripts/README.md`: CLI一覧，実行例，override指定方法
* `docs/phase2/phase2_current.md`: Phase 2の現在地
* `docs/phase2/phase2_tasks.md`: Phase 2の採択判断と根拠
* `docs/phase3/phase3_current.md`: Phase 3の現在地
* `docs/phase3/phase3_tasks.md`: Phase 3の採択判断と根拠
* `docs/phase4/phase4_current.md`: Phase 4の現在地
* `docs/phase4/phase4_tasks.md`: Phase 4の採択判断と根拠
* `docs/codex/codex_workflow.md`: Codex作業，review result，commit／push方針
