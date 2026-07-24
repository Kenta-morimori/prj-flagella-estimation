# Phase 4 Tasks

このファイルは Phase 4 の採択済みタスクを管理する。

チェックボックスは review PASS 後にのみ更新する。

## Phase 4.1: dataset loading contract

### P4-1-001: Issue #146 Phase 3 common clip dataset loader smoke test

- status: complete
- source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/146`
- branch: `feature/phase4-146-clip-loader-smoke`
- goal: Phase 3 common clip dataset (`.npy` clips + `clip_metadata.jsonl`) を Phase 4 が再検出なしに読めることを，学習前の smoke test として確認する。
- result:
  - `src/flagella_estimation/phase4/dataset.py` に Phase 3 common clip dataset loader と audit helper を追加した。
  - loader は `manifest.json`，`clip_metadata.jsonl`，`split_summary.csv`，`.npy` clip を読み，label / split / `group_key` / shape / dtype の contract を検査する。
  - `group_key` leakage，clip frame count mismatch，crop shape mismatch を loader 側で検出できるようにした。
  - training loop / baseline classifier は次Issueへ残した。
- acceptance criteria:
  - [x] Phase 4 が Phase 3 output を再検出なしに読み込める。
  - [x] loader-level test が PASS する。
  - [x] split leakage と label/class count の基本監査が PASS する。
  - [x] baseline classifier へ進むための入力 contract が文書化されている。
- verification:
  - `uv run pytest -q tests/test_phase4_clip_dataset_loader.py tests/test_phase3_gt_passthrough_pipeline.py tests/test_phase3_clip_metadata_schema.py`
- docs:
  - `docs/phase4/phase4_current.md`
  - `docs/phase4/phase4_tasks.md`

## Phase 4.2: diagnostic baseline

### P4-2-001: Issue #148 common clip baseline classifier

- status: complete
- source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/148`
- branch: `feature/phase4-148-baseline-classifier`
- goal: Phase 3 common clip datasetを入力に，既存 grouped splitを維持した最小baseline学習・評価導線を実装する。
- result:
  - 時空間統計特徴量と standardization付き nearest-centroid classifierを `src/flagella_estimation/phase4/` に追加した。
  - `dataset v1`, `n_flagella=1,2,3`, `0.5 s`, `non_overlap`, `qc.status=pass` のfreeze条件を検査する。
  - `model.npz`, `metrics.json`, `predictions.csv`, `confusion_matrix.csv`, `manifest.json`, `run.log` を保存する。
  - 実pilot 18 clips / 9 groupsでCLI smokeを実行し，train / validation / testのartifactを確認した。
  - baseline性能は診断値であり，最終モデル採択や必要独立run数決定には使わない。
- acceptance criteria:
  - [x] Phase 3 common clip datasetからbaselineを学習・評価できる。
  - [x] Phase 3の `group_key` splitを再分割せず維持する。
  - [x] metrics，predictions，confusion matrix，model，manifestを保存する。
  - [x] dataset freeze外のversionを拒否する。
  - [x] library-level testsと実pilot CLI smokeがPASSする。
- verification:
  - `uv run pytest -q tests/test_phase4_baseline_classifier.py tests/test_phase4_clip_dataset_loader.py`
  - `uv run python scripts/04_phase4/train_baseline_classifier.py config=conf/phase4/baseline_v1.yaml dataset_dir=outputs/2026-07-24/001500/phase3_gt_passthrough_v1 output_dir=outputs/2026-07-24/142055/phase4_baseline_v1`
- docs:
  - `docs/phase4/phase4_current.md`
  - `docs/phase4/phase4_tasks.md`

## Phase 4.3: grouped learning curve

### P4-3-001: Issue #150 grouped learning curve evaluator

- status: complete
- source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/150`
- branch: `feature/phase4-150-grouped-learning-curve`
- goal: `track.group_key`を独立sample単位とするclass-balanced learning curve evaluatorを実装し，#129の必要run数判断へartifactを渡す。
- result:
  - 同一groupのclip特徴を平均し，1 groupを1 feature vectorとして扱う。
  - train+val pool内でclassごとにdisjointなtrain / holdout groupsを作り，test splitを最終評価用に保護する。
  - seed固定のbalanced group subsampling，group-level metrics，class recall，empirical percentileを出力する。
  - 全27 v1 candidateで`k=1..6`，各class 2 holdout，100 repeatsを評価した。
  - 現行pseudo v1 / diagnostic baselineでは`k>=4`でmacro F1 meanとempirical intervalが`1.0`になった。
  - 一般的な必要run数の採択は #129 に残し，protected test各class 1 groupという制約を明記した。
- acceptance criteria:
  - [x] x軸がclassごとのunique `group_key` 数である。
  - [x] train / holdout / protected testのgroup leakageを許さない。
  - [x] deterministic seedで再現できる。
  - [x] curve，summary，confusion，manifest，run logを保存する。
  - [x] fixture testsと全27 candidate CLI smokeがPASSする。
- verification:
  - `uv run pytest -q tests/test_phase4_learning_curve.py tests/test_phase4_baseline_classifier.py tests/test_phase4_clip_dataset_loader.py`
  - `uv run python scripts/04_phase4/evaluate_learning_curve.py config=conf/phase4/grouped_learning_curve_v1.yaml dataset_dir=outputs/2026-07-24/143640/phase3_gt_passthrough_v1_full_candidates output_dir=outputs/2026-07-24/144032/phase4_grouped_learning_curve_full_candidates learning_curve.repeats=100`

## Phase 4.4: dataset freeze gate

### P4-4-001: Issue #128 machine-readable dataset freeze audit

- status: complete
- source issue: `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/128`
- branch: `feature/phase4-128-dataset-freeze-audit`
- goal: #128で合意したdataset mixing / versioning規則を，Phase 4 training / learning curve前に実行できるmachine-readable gateへ接続する。
- result:
  - `DatasetFreezePolicy`と全違反を集約するaudit APIを追加した。
  - manifest，metadata，loader auditを使い，version / model / render / torque / class / clip / QC / group provenanceを検査する。
  - Phase 3 manifestの`input_dataset`からPhase 2 `dataset_manifest.json`と選択runの解決済みconfigを辿り，SHA-256と物理regimeを検査する。
  - baseline trainingとgrouped learning curveが同じstrict freeze validatorを通るようにした。
  - PASS / FAILどちらでも`freeze_audit.json`, `manifest.json`, `run.log`を保存するCLIを追加した。
  - 全27 v1 candidate (`54 clips`) のauditがPASSした。
  - 全27 runのsource configを検証し，0 errors / 0 warningsでPASSした。
- acceptance criteria:
  - [x] 主要条件変更の4分類とmixing規則が文書化されている。
  - [x] `model_id`, `render_id`, `dataset_version`, `group_key`のMVP規則を自動検査する。
  - [x] torque / Brownian / RUN-TUMBLE / `n_flagella=4`のMVP除外方針をmachine-readable policyに記録する。
  - [x] training / learning curveがfreeze外datasetを拒否する。
  - [x] PASS / FAIL regression testsと実dataset CLI smokeがPASSする。
- verification:
  - `uv run pytest -q tests/test_phase4_dataset_freeze.py tests/test_phase4_learning_curve.py tests/test_phase4_baseline_classifier.py tests/test_phase4_clip_dataset_loader.py`
  - `uv run python scripts/04_phase4/audit_dataset_freeze.py config=conf/phase4/dataset_freeze_v1.yaml dataset_dir=outputs/2026-07-24/143640/phase3_gt_passthrough_v1_full_candidates output_dir=outputs/2026-07-24/153000/phase4_dataset_freeze_audit_source_verified`
