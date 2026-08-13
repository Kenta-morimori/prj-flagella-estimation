# Phase 2 Issue #61: 2010 torque-linked `dt_star` short-screen campaign

## Scope

この契約は2010 projectの実験専用`1 tau` campaignを定義する。2010 project defaultのlegacy `tau_s=1 s`、#183/ADR 0016のtracking-reference結果、2015 project、0.5秒run、dataset採択を変更しない。

初期screenのconditionはper-flagellum torque `1e-21, 2.5e-20, 1e-19, 1.2e-18 N m`と`dt_star=1e-3, 1e-4`の直積である。`1e-5`は、初期screenで候補を絞った後に実行するformal referenceであり、初期screenの実行対象に含めない。runnerは次を実行前に検証し、違反時はconditionを実行しない。

```text
abs(motor.torque_Nm) = reference_torque_Nm = torque_for_forces_Nm
tau_s = eta * b^3 / T
dt_internal_s = dt_star * tau_s
```

Brownianは無効、seedと2010 projectのmotor/geometryはbase configどおり固定する。

## 実行

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/2010_project_torque_dt_initial_screen.yaml \
  dry_run=true
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/2010_project_torque_dt_initial_screen.yaml
```

実行者はclean worktreeで実行する。出力はJST timestamp付き`outputs/phase2_torque_dt_stability/2010_project/initial_screen/YYYY-MM-DD/HHMMSS/`となる。再集計だけを行う場合は次を使い、simulationを再起動しない。

```bash
uv run python scripts/02_phase2_analysis/analyze_2010_torque_dt_stability.py \
  --run-dir <campaign-root>
```

## First-read QC

`qc_summary.json`、`summary.csv`、`dt_comparison.csv`、`torque_similarity.csv`、`run_manifest.json`をこの順に読む。各conditionはcompact `run_summary.json`、`performance.json`、`state_archive.npz`を持つ。archiveは`t/tau = 0, 0.005, ..., 1.0`の201点でなければ比較不能である。

安全gateは完走、finite、body/non-body shape、motor action-reaction residual（force/torqueとも`<=1e-8`）である。referenceまたはcandidateの安全gate/archiveが失敗したdt比較は`not_comparable`であり、PASSに丸めない。

初期screenでは`1e-4`をscreen comparator、`1e-3`を境界screenとして同一torque内の回転方向、回転量差、姿勢差、COM並進除去bead差を`diagnostic_only`で記録する。初期screenでは数値的一致のPASS/FAILを宣言しない。`1e-4`を暫定候補とする場合に限り、絞り込んだconditionで`1e-5`をformal referenceとして実行し、回転量差`<=10%`、姿勢差`<=5°`、RMS`<=0.1 b`、最大差`<=0.25 b`を正式に判定する。速度、COM変位、姿勢変化は記録するが、1 tauだけでの絶対的遊泳安定性の採否には使わず`diagnostic_review_required`とする。

異torqueの無次元相似性は`dt_star=1e-5`間で記録する診断であり、`1e-4`採否gateではない。初期screenでは`1e-4` screen comparator間の暫定診断として保存するが、formal referenceによる相似性評価と解釈しない。

## Entry-point responsibility

simulation条件の展開・実行は`conf/phase2_multi_run/2010_project_torque_dt_initial_screen.yaml`と標準の`run_multi_run.py`を入口とする。Issue専用の実装は`src/sim_swim/analysis/`に置き、`scripts/01_simulate_swimming/`へ専用runnerを追加しない。実行後のQC再集計・heatmap/replayなどの結果整理は、raw runを再起動しない`02_phase2_analysis`の専用entrypointが担う。

CLI責務は次に固定する。`01_simulate_swimming/`は単発simulation、`run_sweep.py`、`run_multi_run.py`だけを置く。planner・結果解析・plot・既存archiveのreplayは`02_phase2_analysis/`、Phase 2 behavior datasetとPhase 3 clip datasetの作成・replayは`03_dataset_building/`に置く。`scripts/`は薄いentrypointとし、実装本体は`src/`に置く。旧CLI pathには互換wrapperを残さず、active documentationとtestsを新pathへ更新する。
