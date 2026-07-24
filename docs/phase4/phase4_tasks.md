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
