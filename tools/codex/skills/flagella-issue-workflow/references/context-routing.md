# Context Routing

この reference は、issue 単位作業で読む文書を絞るための routing 表である。

## Always Check First

- `git status --short --branch`
- user request の最新内容
- 対象 issue / PR の本文とコメント

## Repository Policy

読む条件:

- branch / commit / push / PR / review_result の扱いを判断する。
- 完了条件や禁止事項を確認する。
- Phase 2 の目視レビュー要否を判断する。

読む文書:

- `AGENTS.md`

注意:

- `AGENTS.md` は正本なので、skill や docs と矛盾した場合は `AGENTS.md` を優先する。

## Project State

読む条件:

- phase の現在地を確認する。
- task が Phase 2.6 / 2.7 / 2.8 など、どの流れに属するか判断する。

読む文書:

- `docs/phase2/phase2_current.md`
- `docs/TASK_MAP.md`
- 必要な場合のみ `docs/PROJECT_PLAN.md`

読み方:

- まず `phase2_current.md` の現在地と参照先だけ読む。
- `PROJECT_PLAN.md` は全体地図や phase-level context が必要な場合だけ、該当 phase の status と直近の progress に絞って読む。
- 全履歴を毎回読み直さない。

## Accepted Tasks

読む条件:

- 既存採択済み task の status、acceptance criteria、branch、source issue を確認する。
- task checkbox を更新してよいか判断する。

読む文書:

- `docs/phase2/phase2_tasks.md`

読み方:

- `rg -n "source issue|P2-|status|acceptance criteria|current diagnostic notes" docs/phase2/phase2_tasks.md`
- 該当 task 周辺だけ読む。

## Phase 2 Detail Docs

読む条件:

- 物理モデル、simulation 条件、既存診断結果、代表出力を確認する。

主な文書:

- `docs/phase2/phase2_6_helix_retention_gate.md`
- `docs/phase2/phase2_6_torque_transmission_model_evaluation.md`
- `docs/phase2/phase2_7_bundling_stability_plan.md`

読み方:

- issue が Phase 2.7 なら Phase 2.7 doc を優先し、Phase 2.6 docs は代表条件や前提が必要な範囲だけ読む。
- 過去の失敗条件を探すときは `outputs/` ではなく、まず docs と `docs/codex-runs/*/review_result.json` を読む。

## Legacy Prompts

読む条件:

- `AGENTS.md` / `docs/` へ移行済みか確認する必要がある。
- 古い acceptance criteria や original prompt の根拠が必要。

読む文書:

- `prompts/00_project_context.md`
- `prompts/01_repo_rules.md`
- relevant `prompts/phase*/`

注意:

- Historical prompts are not source of truth. 現在の正本は `AGENTS.md` と `docs/`。legacy prompt は履歴・移行確認の補助として扱う。

## Codex Run Logs

読む条件:

- 直近 run の PASS / FAIL、commit hash、残 issue、出力先を確認する。
- 既存 PR / branch の作業を引き継ぐ。

探し方:

- `rg -n "source issue|pull_request_url|commit_hash|blocking_issues|next_actions|Issue #<num>|phase2_<task>" docs/codex-runs`

読み方:

- まず `review_result.json` を読む。
- 詳細な reasoning が必要な場合だけ同じ run の `work_log.md` を読む。
- `work_log.md`、長い logs、CSV、生成物はデフォルトで全文表示しない。
