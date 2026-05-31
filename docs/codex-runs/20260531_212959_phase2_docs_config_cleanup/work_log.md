# Phase 2 docs/config cleanup work log

- `docs/phase2/` 配下のユーザー向け文書を `phase2_*` 接頭辞のファイル名へ変更した。
- 旧ファイル名への参照を `docs/PROJECT_PLAN.md`, `docs/TASK_MAP.md`, `docs/planning/phase2_task_proposals.md`, 既存Codexログで更新した。
- `AGENTS.md` に，`docs/phase*/` 配下の Phase 接頭辞命名と，Phase 2 baseline config で default-only `stiffness_scales` を省略する方針を明文化した。
- `conf/sim_swim.yaml` から全値が既定値 `1.0` の `stiffness_scales` ブロックを削除した。
- `uv run pytest` で 130 tests pass を確認した。
