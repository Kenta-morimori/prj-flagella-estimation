# Completion Policy

完了条件と merge 後の引継ぎの正本は `AGENTS.md` と `docs/codex/codex_workflow.md` である。
この reference は、必要な記録先だけを示す。

- run record: `docs/codex-runs/YYYYMMDD_HHMMSS_<phase>_<task-id>/review_result.json`
- required review fields と local PASS: `docs/codex/codex_workflow.md` の「Completion policy」「Review result format」
- `commit → push → source Issue を参照する PR → 初回完了報告`: `AGENTS.md`
- PR マージ後の source Issue 更新と残作業の提示: `docs/codex/codex_workflow.md` の「PRマージ後のIssue引継ぎ」

PR 不要の明示指定または PR 作成失敗時だけ、初回完了報告の順序に例外を設ける。新規 Issue / sub-issue は、
残作業の範囲・受入条件・execution target を提示し、ユーザー承認を得てから作成する。
