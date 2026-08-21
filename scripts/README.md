# scripts

このディレクトリには、リポジトリで直接実行する CLI を置いています。

## 01_simulate_swimming

Phase 2 の 3D 遊泳シミュレーションと条件 sweep を扱います。解析、plot、replay、dataset作成は `03_dataset_building/` の共通CLIに置きます。

### 単発シミュレーション

```bash
uv run python -m scripts.01_simulate_swimming
```

設定はdefaultで `conf/sim_swim_2010.yaml` を読みます。一時的な変更や設定ファイル指定は `KEY=VALUE` 形式で指定します。

```bash
uv run python -m scripts.01_simulate_swimming \
  config=conf/sim_swim_2010.yaml \
  time.duration.value=0.05 \
  time.duration.unit=s \
  motor.torque_Nm=2.0e-20 \
  time.integration.dt_star=1.0e-4
```

`--config`、`--duration-s`、`--fps-out`、`time.duration_s`、`time.dt_star` は既存互換用に残しています。新規コマンド例では `config=...`、`time.duration.value=...`、`time.duration.unit=...`、`time.integration.dt_star=...`、`output_sampling.fps_out_2d=...` を使います。

時間正規化の契約:

- `time.duration.unit` は `s` または `tau`。
- 2010 project profile は legacy 互換のため `tau_s=1.0` 固定。
- 2010 paper / 2015 profiles は `tau_s = viscosity_Pa_s * b_m^3 / abs(motor.reference_torque_Nm)`。
- `time.dt_s` は deprecated な出力・記録間隔で、内部積分刻みではありません。
- `total_steps = ceil(duration_tau / dt_star)`。`step_summary.csv` はstep開始時刻を記録するため最終state時刻は指定durationをわずかに超える場合があります。

### Issue #190: 2010 project torque連動time-scaleの短時間screen

既存の2010 project defaultは`tau_s=1 s`固定のまま保持する。Issue #190のplanだけは
`time.scale_policy=reference_torque`により、各conditionで`motor.torque_Nm`、
`reference_torque_Nm`、物性scale torqueを同一値にし、`tau_s=eta*b^3/T`を使う。
これは短時間の数値安定性実験であり、既存baseline、dataset、2015 profileを採択しない。

```bash
uv run python scripts/03_dataset_building/analyze_dataset.py --analysis-kind plan-2010-torque-dt \
  --config conf/phase2_sweeps/2010_project_torque_linked_dt_stability.yaml \
  --output /private/tmp/issue190_campaign_plan.json
```

初期screenは`1 tau`であり、4 torque × 3 `dt_star`の12条件をsimulationなしで展開する。
`dt_star=1e-4`が候補default、`1e-3`が境界screen、`1e-5`が比較referenceである。

Issue #61の初期screenは`dt_star=1e-3,1e-4`だけを標準multi-runで実行する。
`1e-5`は初期screenで絞った候補を正式に検証するための後続referenceであり、初期screenの
`diagnostic_only`結果を正式採否に読み替えない。実行と既存campaignの再集計は次です。

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/2010_project_torque_dt_initial_screen.yaml \
  dry_run=true
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/2010_project_torque_dt_initial_screen.yaml
uv run python scripts/03_dataset_building/analyze_dataset.py --analysis-kind 2010-torque-dt \
  --run-dir <campaign-root>
```

campaign root の `qc_summary.json` を最初に読みます。`summary.csv` はcondition safety、
`dt_comparison.csv` は同一torque内の`1e-3` / `1e-4` screen一致度、`torque_similarity.csv` は
異torque無次元相似性の暫定診断です。`1e-5` formal reference比較は、初期screenで候補を絞った後に別campaignで行います。1 tauの遊泳速度・姿勢は記録し、単独では採否に使いません。

既存campaignの形状・後方軸整列・束化診断をtorque × `dt_star` heatmapへまとめる場合は、simulationを再起動せず次を使います。

```bash
uv run python scripts/03_dataset_building/analyze_dataset.py --analysis-kind visualize-2010-torque-dt \
  --run-dir <campaign-root>
```

`analysis/torque_dt_visuals/torque_dt_feature_heatmaps.png`と、同じ値を持つCSVを出力します。最初のpanelは誤差の合計ではなく、完走・有限性・形状・motor action-reactionを含む必須QCのAND（`PASS` / `FAIL`）です。ほかのpanelは、flag bond/hook長の最大相対誤差、べん毛軸の最大平均ずれ、後方軸への最小投影、最終frameの束参加率、束半径の最大値を示します。`1 tau`のheatmapは短時間診断であり、論文モデルの長時間束化や`dt_star`の正式採用を表すものではありません。

同じ8条件を同一の無次元時刻で3D replay gridへ並べる場合は次を使います。

```bash
uv run python scripts/03_dataset_building/replay_dataset.py \
  --run-dir <campaign-root> \
  --mode render-only \
  --target-frame-count 101 \
  --max-panels-per-grid 8
```

すべての3D renderは、policyにかかわらず` t = … τ (… s, … steps)`を表示する。motor torque・scale・QCなどの条件情報はこの時刻行の後へ続けて表示する。

主な出力は `outputs/YYYY-MM-DD/HHMMSS/` 配下に作成され、`manifest.json` と `run.log` に実行条件が記録されます。Issue #61 torque-dt screenのrun rootは`outputs/YYYY-MM-DD/HHMMSS/phase2_issue61_2010_project_torque_dt_initial_screen/`です。terminalと`run.log`へ全8条件の開始・完了/失敗と経過時間が出力され、実行中のconditionは`run_manifest.json`の`execution.status`で確認できます。

Issue #193は、2010 torque-linked条件で同じ物理`0.5 s`を回すcostだけを測ります。無次元時間の同等性や`dt_star`採択は判断しません。

Issue #199は、τ固定controlとtorque-linked τを同じ物理`0.05 s`で比較する12条件screenです。実行後は同じcampaign rootを再解析・可視化・replayへ渡します。

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/2010_project_tau_policy_torque_dt_0p05s.yaml \
  dry_run=true
uv run python scripts/03_dataset_building/analyze_dataset.py --analysis-kind 2010-torque-dt \
  --config conf/phase2_multi_run/2010_project_tau_policy_torque_dt_0p05s.yaml \
  --run-dir <campaign-root>
uv run python scripts/03_dataset_building/analyze_dataset.py --analysis-kind visualize-2010-torque-dt \
  --run-dir <campaign-root>
```

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/2010_project_torque_fixed_real_time_0p5s_performance.yaml
uv run python scripts/03_dataset_building/analyze_dataset.py --analysis-kind 2010-fixed-performance \
  --run-dir <campaign-root>
```

`performance_summary.csv`には、各torqueの`duration_tau`、総step数、simulation-loop wall time、steps/sと必須safety QCを出力します。

### Sweep

開発用・診断用の task-specific sweep は `run_sweep.py` を使います。対象 profile は `conf/phase2_sweeps/` に置きます。

```bash
uv run python scripts/01_simulate_swimming/run_sweep.py \
  config=conf/phase2_sweeps/shape_stability_grid.yaml \
  dry_run=true \
  sample_limit=3
```

実行時は `[1/3] shape_stability_grid ...` のように進捗が標準出力へ表示されます。

```bash
uv run python scripts/01_simulate_swimming/run_sweep.py \
  config=conf/phase2_sweeps/shape_stability_grid.yaml \
  mode=first-second-grid \
  time.duration.value=0.001 \
  time.duration.unit=s \
  motor.torque_Nm=0 \
  first_second_spring_scales=1 \
  output_dir=/private/tmp/phase2_smoke \
  overwrite=true
```

profile の既定値は `KEY=VALUE` で上書きできます。sweep の標準 summary は出力先の `summary.csv` です。`--config` などの legacy option 形式も互換用に残しています。

### Multi-Run

ユーザが複数条件を一度で実行し、その結果から replay / plot へ進む入口は `run_multi_run.py` です。profile は `conf/phase2_multi_run/` に置き、1つの config を run / plot / replay で共用します。campaign 単位で `run.log`、`manifest.json`、`run_manifest.json`、`summary.csv` を残します。

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/latest_model_torque_shape_stability.yaml \
  dry_run=true
```

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/latest_model_torque_shape_stability.yaml \
  sweep.axes.torque.values=[1.5e-20,2.0e-20,2.5e-20] \
  time.duration.value=1.0 \
  time.duration.unit=s \
  overwrite=true
```

`output.timestamp_subdir=false` の multi-run profile では、`output.base_dir` がそのまま run root になります。再実行で同じ run root を置き換える場合だけ `overwrite=true` を明示します。

主な profile:

| profile | 用途 |
| --- | --- |
| `conf/phase2_sweeps/motor_scale.yaml` | motor-local scale の sweep |
| `conf/phase2_sweeps/single_flagellum_torque.yaml` | single flagellum torque 条件評価 |
| `conf/phase2_sweeps/bundling_alignment.yaml` | 複数べん毛の helix axis alignment 評価 |
| `conf/phase2_sweeps/shape_stability_grid.yaml` | hook / proximal flagellum を含む shape stability grid |
| `conf/phase2_sweeps/torque_distribution_grid.yaml` | #97 用 torque distribution 2x2 比較 |
| `conf/phase2_sweeps/hook_overstretch.yaml` | 旧名互換 profile |
| `conf/phase2_multi_run/spring_formulation_motor_off.yaml` | Issue #163 motor-off `legacy` / `fene_fraenkel` 比較 |
| `conf/phase2_multi_run/spring_formulation_motor_on.yaml` | Issue #163 motor-on RUN `legacy` / `fene_fraenkel` 比較 |

新規の user-facing 実行では `shape_stability_grid.yaml` を正本として使います。
`hook_overstretch.yaml` は historical alias であり，既存メモや過去コマンドの互換用です。

利用可能な sweep profile を CLI から確認する場合:

```bash
uv run python scripts/01_simulate_swimming/run_sweep.py \
  list_canonical_profiles=true
```

個別 profile の `role` / `canonical` / 推奨 heatmap を確認する場合:

```bash
uv run python scripts/01_simulate_swimming/run_sweep.py \
  config=conf/phase2_sweeps/shape_stability_grid.yaml \
  describe_profile=true
```

### 2010 spring formulation比較

`potentials.spring.formulation` は `legacy` と `fene_fraenkel` を受け付けます。selectorを省略した既存configは過去互換の `legacy` として読み込みます。単発で論文準拠式を選ぶ場合:

```bash
uv run python -m scripts.01_simulate_swimming \
  config=conf/sim_swim_2010.yaml \
  potentials.spring.formulation=fene_fraenkel \
  time.duration.value=0.01 \
  time.duration.unit=s \
  time.integration.dt_star=1.0e-4
```

Issue #163 のdefault判定用長時間比較は、motor-offとmotor-on RUNを同じseedで実行します。

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/spring_formulation_motor_off.yaml

uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/spring_formulation_motor_on.yaml
```

各コマンドが出力したrun rootの `summary.csv` を解析CLIへ渡します。

```bash
uv run python scripts/03_dataset_building/analyze_dataset.py --analysis-kind spring-formulations \
  motor_off_summary=<motor-off-run-root>/summary.csv \
  motor_on_summary=<motor-on-run-root>/summary.csv
```

解析runには `force_extension.csv` / `force_extension.png` と、自動採否の `default_decision.json` / `default_decision.md` が出力されます。採否は両runの完走、有限値、body/nonbody strict shape gateを使い、条件を満たす場合だけ `fene_fraenkel` を選びます。この比較に動画目視は要求しません。

### 2015 refined Stage A

Issue #168はmotor-off pilotを先に実行し、閾値を固定してからmotor-onへ進みます。

```bash
uv run python scripts/01_simulate_swimming/run_sweep.py \
  config=conf/phase2_sweeps/2015_stage_a_motor_off.yaml

uv run python scripts/03_dataset_building/analyze_dataset.py --analysis-kind 2015-stage-a \
  --motor-off-run <motor-off-run-root>
```

`dt_star=1e-4`参考比較、motor-on command、brace分岐、replay、目視項目は
`docs/phase2/phase2_168_2015_stage_a_validation.md`を参照してください。motor-off pilotから具体的閾値を
lockする前にmotor-on採否は行いません。

Issue #199 Task D の2015 project τ-linked torque × `dt_star` screen は、reference/candidateの
2 runを完了後に統合します。実行条件と診断限定の扱いは
`docs/phase2/phase2_199_tau_policy_torque_dt_campaign_contract.md` を正本とします。

### Heatmap

sweep summary から heatmap を作る場合は `plot_heatmap.py` を使います。

```bash
uv run python scripts/03_dataset_building/analyze_dataset.py \
  config=conf/phase2_sweeps/shape_stability_heatmap.yaml \
  summary_csv=/private/tmp/phase2_smoke/summary.csv \
  mode=first-second-grid
```

heatmap profile は出力先を固定しません。`output_dir` を省略すると、`summary_csv` と同じ directory の `plots/` へ出力します。明示した場合はその出力先を使います。
`shape_stability_heatmap.yaml` は `mode=position-only-grid` などの実行時 override で対象 grid を切り替えられます。

generic multi-run の summary plot も `plot_heatmap.py` から行います。同じ config をそのまま使います。`plot.default_y_axis` が未設定の profile では heatmap ではなく 1 軸 line plot を出します。

```bash
uv run python scripts/03_dataset_building/analyze_dataset.py \
  config=conf/phase2_multi_run/latest_model_torque_shape_stability.yaml
```

`output.timestamp_subdir=false` の multi-run profile では、`summary_csv` / `run_dir` を省略すると `output.base_dir/summary.csv` を読み、`output.base_dir/plots/` へ出力します。

主な profile:

| profile | 用途 |
| --- | --- |
| `conf/phase2_sweeps/motor_scale_heatmap.yaml` | motor scale sweep の collapse heatmap |
| `conf/phase2_sweeps/dt_star_torque_heatmap.yaml` | `dt_star` x torque heatmap |
| `conf/phase2_sweeps/local_scale_mode_heatmap.yaml` | local scale mode x torque heatmap |
| `conf/phase2_sweeps/shape_stability_heatmap.yaml` | shape stability grid heatmap |
| `conf/phase2_multi_run/latest_model_torque_shape_stability.yaml` | generic multi-run summary plot / replay metadata |
| `conf/phase2_sweeps/hook_overstretch_heatmap.yaml` | 旧名互換 heatmap profile |

heatmap も `shape_stability_heatmap.yaml` を正本として使います。
`hook_overstretch_heatmap.yaml` は historical alias です。

利用可能な heatmap profile を CLI から確認する場合:

```bash
uv run python scripts/03_dataset_building/analyze_dataset.py \
  list_canonical_profiles=true
```

### Replay Render

既存 sweep 出力の `summary.csv`、`run_manifest.json`、各 condition directory の `state_archive.npz` から、再シミュレーションなしで比較 plot / 2D・3D replay を生成する場合は `03_dataset_building/replay_dataset.py` を使います。generic multi-run 出力でも同じ CLI を使います。

```bash
uv run python scripts/03_dataset_building/replay_dataset.py \
  config=conf/phase2_multi_run/latest_model_torque_shape_stability.yaml \
  overwrite=true
```

`--mode plot-only` は metrics CSV / PNG のみ、`--mode render-only` は 3D grid movie のみ、`--mode both` は両方を生成します。
`output.timestamp_subdir=false` の multi-run profile では、`run_dir` / `input_dir` を省略すると `output.base_dir` を読み、`output.base_dir/replay/` へ出力します。legacy 互換として `summary_csv=...` や `--input-dir ... --output-dir ...` も引き続き使えます。

## 共通 dataset workflow

`scripts/02_phase2_analysis/` は PR #201 で廃止した。既存 raw run / dataset の再利用は次の4操作だけを使う。

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py config=<campaign.yaml>
uv run python scripts/03_dataset_building/build_dataset.py config=<campaign.yaml> run_dir=<run-root>
uv run python scripts/03_dataset_building/replay_dataset.py --dataset-dir <dataset-dir> --view 3d+2d
uv run python scripts/03_dataset_building/analyze_dataset.py --dataset-dir <dataset-dir> --view both
```

`inspect_run.py` は `run_summary.json` を入口にした限定診断用であり、large `step_summary.csv` を
全件読みする代替ではない。

### CLI導線・lifecycle一覧

| 分類 | 公開command | 用途 / owner | 実装 | lifecycle |
| --- | --- | --- | --- | --- |
| canonical | `01_simulate_swimming.py` / `run_sweep.py` / `run_multi_run.py` | simulation・sweep・multi-run | `src/sim_swim/sim/`, `src/sim_swim/analysis/sweeps/` | canonical |
| canonical | `build_dataset.py` / `build_clip_dataset.py` | behavior / clip dataset作成 | `src/sim_swim/analysis/`, `src/flagella_estimation/phase3/` | canonical |
| canonical | `analyze_dataset.py --dataset-dir ...` | datasetの3D分布・2D識別性分析 | `src/sim_swim/analysis/behavior_dataset_*`, `src/flagella_estimation/phase3/feature_comparison.py` | canonical |
| canonical | `replay_dataset.py` | raw run / datasetの非再simulation replay | `src/sim_swim/analysis/phase2_replay.py`, `behavior_dataset_replay.py` | canonical |
| canonical | `analyze_dataset.py --analysis-kind heatmap` / `run-summary` | profile heatmap / compact run summary | `src/sim_swim/analysis/common_analysis.py`, `run_summary.py` | canonical |
| active diagnostic | `analyze_dataset.py --analysis-kind 2010-torque-dt`, `2010-fixed-performance`, `visualize-2010-torque-dt`, `partial-generic` | #61 / #193 / #199 torque・`dt_star` campaignとpartial集計 | `torque_dt_stability_campaign.py`, `torque_dt_stability_visuals.py`, `partial_generic_multi_run.py` | active diagnostic |
| active diagnostic | `analyze_dataset.py --analysis-kind plan-2010-torque-dt`, `plan-reference-torque` | #61 / #183 campaign plan | `torque_dt_stability.py`, `reference_torque_comparison.py` | active diagnostic |
| active diagnostic | `analyze_dataset.py --analysis-kind diagnose-158` / `inspect_run.py` | #158 failure・限定run診断 | `phase2_158_diagnostics.py`, `run_summary.py` | active diagnostic |
| active diagnostic | `analyze_dataset.py --analysis-kind task-d-2015`, `phase-seed-2d` / `evaluate_clip_dataset_features.py` | #199 Task D・v1 r2 feature確認 | `task_d_2015_tau_linked.py`, `phase_seed_2d_replay.py`, `src/flagella_estimation/phase3/feature_comparison.py` | active diagnostic |
| result reproduction | `analyze_dataset.py --analysis-kind spring-formulations`, `2015-stage-a` | #163 / #168の既存結果・診断再現 | `spring_formulations_workflow.py`, `stage_a_2015_workflow.py` | historical result reproduction |
| result reproduction | `render_sample.py` / `replay_clip_dataset.py` | Phase 3 dataset・clipの再現 | `src/flagella_estimation/phase3/` | historical result reproduction |

Issue専用の分析kindは、campaignのlive contractから実行する。追加・削除時は
`docs/codex/issue_script_lifecycle.md` の `promoted` / `deleted` /
`retained-temporarily` 分類を更新し、`tools/codex/check_issue_script_lifecycle.py` を確認する。

以下の command 例は、結果の来歴を説明するためだけに残す。

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/flagella_count_behavior_v0.yaml \
  dry_run=true sample_limit=3
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/flagella_count_behavior_v0.yaml
uv run python scripts/03_dataset_building/build_dataset.py \
  config=conf/phase2_multi_run/flagella_count_behavior_v0.yaml
```

標準設定:

| config | 用途 |
| --- | --- |
| `conf/phase2_multi_run/flagella_count_behavior_v0.yaml` | Issue #71 の診断用 dataset v0。run / heatmap / replay / dataset 作成で共通に使う |

Issue #71 の診断用 dataset v0 は次で実行します。36 sample の長時間 run なので、まず `dry_run=true` や `sample_limit=1` で確認してください。

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/flagella_count_behavior_v0.yaml
uv run python scripts/03_dataset_building/analyze_dataset.py \
  config=conf/phase2_multi_run/flagella_count_behavior_v0.yaml
uv run python scripts/03_dataset_building/replay_dataset.py \
  config=conf/phase2_multi_run/flagella_count_behavior_v0.yaml
uv run python scripts/03_dataset_building/build_dataset.py \
  config=conf/phase2_multi_run/flagella_count_behavior_v0.yaml
uv run python scripts/03_dataset_building/analyze_dataset.py --analysis-kind distributions \
  --dataset-id v0
```

raw condition は `step_summary.csv`、`trajectory.csv`、`state_archive.npz` を保存します。dataset directory から 3D / 2D render を作る場合は次を使います。

```bash
uv run python scripts/03_dataset_building/render_sample.py \
  --dataset-dir outputs/phase2_analysis/flagella_count_behavior/datasets/v0
```

dataset の分布 plot は次で生成します。

```bash
uv run python scripts/03_dataset_building/analyze_dataset.py --analysis-kind distributions \
  --dataset-id v0
```

dataset v1 の `n_flagella=1,2,3` について，2D投影後の運動特徴量と grouped baseline を確認する場合は次を使います。

```bash
uv run python scripts/03_dataset_building/analyze_dataset.py \
  --dataset-dir outputs/phase2_analysis/flagella_count_behavior/datasets/v1 \
  --output-dir outputs/phase2_analysis/flagella_count_behavior/datasets/v1/analysis/projection_2d \
  --overwrite
```

## 03_dataset_building

Phase 2のbehavior datasetとPhase 3のclip datasetを、既存simulation/archiveから独立して作成します。解析とdataset作成は `03_dataset_building/`、Phase 4の学習・評価は `04_phase4/` が担当します。

### v1 r1 3秒runから0.5秒clip datasetを作成

Issue #159 の Phase 3 common clip dataset は、既存の3秒 state archive から生成します。長時間の Phase 2 simulation は再実行しません。canonical `.npy` clip は、べん毛を描画せず、菌体のみを `body_capsule_orthographic_v1` のtracking cropで描画します。菌体は全長2.0 µm、幅1.0 µmの剛体capsuleとして扱い、固定カメラへの並行投影により斜め姿勢では見かけの長さを短縮し、z方向への正対時は円形にします。菌体変形はML入力へ反映せず、QC確認対象として扱います。blur / noise / defocusは実動画条件が未確定のためcanonical v1には含めません。Phase 3 configでは、出力fpsは `output_sampling.fps_out`、画像条件は `render.image_size_px` / `render.pixel_size_um` を正本とします。

画像サイズ、pixel scale、body寸法、intensity等をoverrideした場合は描画条件から決定的なvariant `render_id`を生成します。canonical freeze auditはvariant IDを拒否します。MP4 replayはdataset manifestの描画条件を復元するため、保存済み`.npy`と同じsilhouetteを表示します。

probe は1 classあたり1 runだけ処理します。probe 出力は実行ごとの timestamp path に置き、最終 dataset v1 には使いません。

```bash
uv run python scripts/03_dataset_building/build_clip_dataset.py \
  config=conf/phase3/gt_passthrough_v1_r1_duration_3s_clips.yaml \
  filters.max_per_class=1
```

final build は27 runすべてを処理し、Phase 3 common clip dataset v1 の canonical path として `outputs/phase3_common_clip/datasets/v1/` に出力します。期待値は5 clips/run、合計135 clips、独立group 27です。

```bash
uv run python scripts/03_dataset_building/build_clip_dataset.py \
  config=conf/phase3/gt_passthrough_v1_r1_duration_3s_clips.yaml \
  output_dir=outputs/phase3_common_clip/datasets/v1
```

`outputs/YYYY-MM-DD/HHMMSS/phase3_v1_r1_clip_dataset/` はstaging / probe用途の実行ログ付き出力です。最終採択済みの Phase 4 入力 dataset は `outputs/phase3_common_clip/datasets/v1/` を参照します。`data/` は外部入力や手元データ置き場として残し、今回の生成済み共通clip dataset正本には使いません。

QC確認対象は `dataset_summary.csv`、`qc_summary.csv`、`split_summary.csv`、`clip_metadata.jsonl` です。`n_flagella=3` の first-fail run では、first-failを含むwindowとそれ以降のwindowが `qc_label=diagnostic` / `training_candidate=false` になります。ここでの `diagnostic` は、Phase 2 shape QC の `first_fail_t_s` を含む、またはそれ以後のためtrainingから除外した診断用clipという意味であり、全frameで形状崩壊が目視できることを意味しません。early clipは削除せず、Phase 4側の `freeze.warmup_s=0.0/0.5/1.0` で選択します。

Phase 4 baselineとgrouped learning curveは、どちらも`training_candidate=true`かつ`clip.t_start_s >= freeze.warmup_s`の同じ選択集合だけを学習・集約へ渡します。diagnostic clipはartifactとして残りますが、group featureには混入しません。
v1 r1 では source summary の `use_for_ml_candidate=false` run も除外せず、clip artifact を作成したうえで window QC により diagnostic-only として扱います。

MP4 grid replay は source state archive から同じclip windowを再構成し、defaultでは各panelに3D source replayと2D body-only renderを横並びで表示します。3D側はPhase 2共通3D rendererを使い、べん毛を描画します。2D側はcanonical ML入力と同じ菌体のみのcapsule renderです。titleには `n_flagella`、`run_id`、`clip_index`、time band、QC label が入り、3D status textには first-fail時刻と before/after first-fail status が入ります。

`--max-clips` を省略すると、filterに一致する全clipを処理します。1本のMP4に描画する上限はdefault 12 clipで、対象が12 clipを超える場合は `3d_2d_grid_001.mp4`、`3d_2d_grid_002.mp4` のように分割します。この上限は `--clips-per-video` で変更できます。`--max-clips` は全体の処理件数を明示的に制限したいprobe用途に使います。

```bash
uv run python scripts/03_dataset_building/replay_clip_dataset.py \
  outputs/phase3_common_clip/datasets/v1
uv run python scripts/03_dataset_building/replay_clip_dataset.py \
  outputs/phase3_common_clip/datasets/v1 \
  --qc-label diagnostic
uv run python scripts/03_dataset_building/replay_clip_dataset.py \
  outputs/phase3_common_clip/datasets/v1 \
  --training-candidate true \
  --n-flagella 3 \
  --max-clips 12
```

MP4 panelは必要に応じて選択できます。defaultは `--panel-layout 3d+2d` です。

```bash
uv run python scripts/03_dataset_building/replay_clip_dataset.py \
  outputs/phase3_common_clip/datasets/v1 \
  --panel-layout 3d
uv run python scripts/03_dataset_building/replay_clip_dataset.py \
  outputs/phase3_common_clip/datasets/v1 \
  --panel-layout 2d
```

従来の `.npy` tensor contact sheetを確認したい場合は `--mode contact-sheet` を付けます。

```bash
uv run python scripts/03_dataset_building/replay_clip_dataset.py \
  outputs/phase3_common_clip/datasets/v1 \
  --mode contact-sheet \
  --max-clips 12
```

定性確認点は、early/middle/late clipでtracking cropが破綻していないこと、2D側にflagellaが描画されていないこと、菌体がcapsule状で不自然に伸縮しないこと、3Dと2Dの姿勢・時刻対応が破綻していないこと、pre-first-failとdiagnosticの境界表示がmetadataと一致すること、class/run/time bandが期待どおり見分けられることです。

Phase 4 freeze audit は同じ Phase 3 dataset を直接読みます。`warmup_s` は比較条件として切り替えます。

```bash
uv run python scripts/04_phase4/audit_dataset_freeze.py \
  config=conf/phase4/dataset_freeze_v1_r1.yaml \
  dataset_dir=outputs/phase3_common_clip/datasets/v1
uv run python scripts/04_phase4/audit_dataset_freeze.py \
  config=conf/phase4/dataset_freeze_v1_r1.yaml \
  dataset_dir=outputs/phase3_common_clip/datasets/v1 \
  freeze.warmup_s=0.5
uv run python scripts/04_phase4/audit_dataset_freeze.py \
  config=conf/phase4/dataset_freeze_v1_r1.yaml \
  dataset_dir=outputs/phase3_common_clip/datasets/v1 \
  freeze.warmup_s=1.0
```

## 04_phase4

### Duration / Seed Study

Issue #129 のサンプル時間長・run数検討では、3秒raw run、dataset revision、duration / seed解析を順番に実行します。条件だけを確認する場合は `probe`、full factorial 27 runを実行する場合は `full` を指定します。

```bash
./scripts/04_phase4/run_duration_seed_study.sh probe
caffeinate -i ./scripts/04_phase4/run_duration_seed_study.sh full
```

raw run、dataset revision、duration / seed解析は、実行ごとのJST timestampを持つ同じcampaign rootの `dataset/`、`analysis/` 配下へ保存されます。

```text
outputs/YYYY-MM-DD/HHMMSS/phase4_duration_seed_study/
├── dataset/v1_r1_duration_3s/
├── analysis/phase4_duration_seed_study/
└── replay/
```

同一特徴量のtime plotは `0.25 / 0.5 / 1.0 s` で共通y軸を使います。seed heatmapはdurationを行、`n_flagella` を列とする1枚にまとめ、durationごとに独立したcolor scaleを使います。`x` と `!` はそれぞれQC fail windowとQC fail runを示します。

3秒raw runがすでに存在する場合は、`RUN_DIR`にcampaign rootを指定し、シミュレーションを再実行せずdataset revisionと解析結果だけを生成できます。以下は既存のcanonical campaignを再解析する例です。

```bash
RUN_DIR=outputs/phase2_multi_run/flagella_count_duration_3s_r1
uv run python scripts/03_dataset_building/build_dataset.py \
  config=conf/phase2_multi_run/flagella_count_duration_3s_r1.yaml \
  run_dir="${RUN_DIR}" \
  dataset.output_dir="${RUN_DIR}/dataset/v1_r1_duration_3s"
uv run python scripts/04_phase4/evaluate_duration_seed_study.py \
  config=conf/phase4/duration_seed_study_v1_r1.yaml \
  dataset_dir="${RUN_DIR}/dataset/v1_r1_duration_3s" \
  output_dir="${RUN_DIR}/analysis/phase4_duration_seed_study"
```

3秒raw runを再シミュレーションせず3D比較renderする場合:

```bash
uv run python scripts/03_dataset_building/replay_dataset.py \
  config=conf/phase2_multi_run/flagella_count_duration_3s_r1.yaml \
  run_dir="${RUN_DIR}" \
  mode=render-only \
  max_panels_per_grid=9 \
  overwrite=true
```

27条件は最大9条件ずつ複数pageへ出力されます。各パネルにはconditionとともに `PASS`、または `FAIL <category>@<first_fail_t_s>` が表示されます。
