# Completion Policy

この reference は、issue 対応の完了処理を軽く確実に行うための要約である。正本は `AGENTS.md`。

## Run ID

形式:

`YYYYMMDD_HHMMSS_<phase>_<task-id>`

例:

`20260615_233000_codex_issue50_skill`

保存先:

`docs/codex-runs/<run-id>/`

## Required Files

原則:

- `docs/codex-runs/<run-id>/review_result.json`
- 必要なら `docs/codex-runs/<run-id>/work_log.md`

`review_result.json` は完了状態の正本である。

## PASS Conditions

PASS と言うには、次を満たす。

- requested implementation / documentation change が完了している。
- relevant tests/checks が PASS、または未実行理由が明確。
- review step が完了している。
- `review_result.json` の `status` が `"PASS"`。
- work log / review result が保存されている。
- final state が commit 済み。
- remote access があれば push 済み。
- pushed feature branch なら PR 作成済み。

## FAIL / Diagnostic Conditions

Phase 2 では、FAIL でも有用な診断進捗を commit / push してよい。

例:

- collapse mode を再現した。
- fly-away 条件を保存した。
- hook drift や no_bundle_drive を切り分けた。
- failing test で次の target behavior を定義した。
- 物理モデル差分や数値上の不一致を記録した。

FAIL commit は完了扱いにしない。commit type は `diagnostic`, `wip`, `docs`, `test` 相当の内容にする。

## Review Result Fields

最低限、以下を埋める。

- `status`
- `summary`
- `blocking_issues`
- `non_blocking_issues`
- `tests_reviewed`
- `user_review_required`
- `user_review_command`
- `user_review_outputs`
- `user_review_points`
- `adr_required`
- `adr_reason`
- `commit_type`
- `commit_hash`
- `push_status`
- `pull_request_url`
- `next_actions`

## Commit / Push / PR

- commit message は `type(scope): summary`。
- `main` / `master` へ直接 commit しない。
- feature branch を push したら PR を作る。
- GitHub issue は明示依頼がない限り作らない。
- PR は作ってよいが、merge はしない。
