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
