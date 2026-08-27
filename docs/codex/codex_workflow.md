# Codex Workflow

この文書は，毎回読む必要のない Codex 運用詳細をまとめる。

通常は `AGENTS.md` と対象 task の current doc だけを読み，完了条件・review_result・commit/push/PR・ADR・Cloud review の判断が必要なときだけこの文書を読む。

## Source of truth

`docs/codex-runs/<run-id>/review_result.json` は Codex task 完了状態の正本である。

Task checkboxes，commit message，PR本文，Codex final response は二次記録であり，`review_result.json` と矛盾してはいけない。

`review_result.json` の `PASS` は，ローカル実装・文書・セルフチェックが完了していることを表す。PR作成後の CI と trusted Cloud review（CodexまたはGitHub Copilot）は merge gate であり，PR checklist と GitHub checks で確認する。trusted review 未実施だけを理由に，ローカル完了済みの `review_result.json` を `FAIL` に戻さない。

## Issue execution target

新規IssueはGitHub Issue Formでheavy/runtime execution targetを必須選択する。コード実装と
短時間unit testは通常Macで行い、このtargetは長時間simulation、sweep、render等の実行先を示す。

| Form value | Label | Permitted execution |
| --- | --- | --- |
| `mac_only` | `execution:mac` | Mac local runtimeのみ |
| `cs10_user_run` | `execution:cs10` | Mac実装・短時間check後、Userがcs10でheavy jobを実行 |
| `no_runtime` | `execution:none` | docs / review / workflowのみ |
| `triage_required` | `execution:triage` | read-only triageのみ |

独立conditionが8以上、またはMac見積りwall timeが30分超なら`cs10_user_run`を選ぶ必須候補とする。
Issue作成・編集時のworkflowが`execution:*` labelを同期する。本文のtargetとlabelが不在・不一致、
または`execution:triage`なら、Codexは実装・test・runtimeを開始せず、targetのtriageを依頼する。
既存Issueは一括推測せず、着手時にこのForm項目を追記してtriageする。

## Issue Roadmap metadata

新規IssueはIssue Formの必須`Roadmap category (Milestone)`を選択する。`issue-roadmap-sync` workflowは
Project #8へIssueを登録し、対応Milestoneと作成日（JST）の`Start date`を同期する。
`Planned target date`は任意の`YYYY-MM-DD`入力であり、指定時は`Target date`へ同期する。Issue close時に
Target dateが未設定なら、そのIssueの終了日（JST）を補完する。既に予定日があれば上書きしない。

このworkflowはuser-owned Projectを更新するため、classic PAT（`repo` + `project` scope）を
`PROJECT_AUTOMATION_TOKEN`としてrepository secretへ登録する。tokenが無い場合は明示的に失敗する。
Issue Form外で作成されたIssueや不正な日付は`roadmap:triage`、reopenされたIssueは
`roadmap:needs-review`で明示し、metadata修正まで実装に着手しない。

PR URL，最終PR head SHA，push後の状態など，PR作成後にしか確定しない動的情報を tracked `review_result.json` へ後追い同期するためだけの commit は作らない。これらはPR本文，GitHub checks，最終ユーザー報告に記録する。

## Run ID

形式:

`YYYYMMDD_HHMMSS_<phase>_<task-id>`

例:

`docs/codex-runs/20260530_142233_phase2_0037/review_result.json`

## Completion policy

PASS 完了には以下が必要である。

1. Requested implementation / documentation change が完了している。
2. Relevant tests/checks が PASS，または未実行理由が明確である。
3. Review step が完了している。
4. `docs/codex-runs/<run-id>/review_result.json` が `"status": "PASS"` である。
5. Work log / review result が保存されている。
6. Final state が commit 済みである。
7. Remote access があれば push 済みである。
8. Pushed feature branch なら PR が作成済みである。

FAIL result は完了ではない。ただし Phase 2 では，有用な診断進捗を `diagnostic`, `wip`, `docs`, `test` 相当の commit として保存してよい。

有用な FAIL 例:

* collapse / fly-away / hook drift / no_bundle 条件を再現した。
* failing test で次の target behavior を定義した。
* 物理モデル差分や数値上の不一致を記録した。
* 部分実装で原因範囲を狭めた。

## Test policy

pre-commit hook の既定は lightweight checks とする。

既定コマンド:

* `uv run ruff format --check .`
* `uv run ruff check .`
* `uv run pytest -q -m light`

`light` は commit 時に固定実行しても負担が小さい，短時間・deterministic・library-level の test を指す。初期運用では対象を広げすぎず，明らかに軽い test だけを明示的に marker 付与する。

full pytest は削除しない。以下では `uv run pytest -q` を実行する。

* 物理モデル，geometry，hook，flagella，body，torque，force，potential，hydrodynamics の変更時
* simulation core の変更時
* output format，manifest，CSV schema の変更時
* dataset 生成仕様の変更時
* PR 作成前または merge 前
* GitHub Actions CI

docs-only，planning-only，workflow-only 変更では full pytest を既定要求しない。ただし未実行の場合も，必要なら `review_result.json` に理由を書く。

hook で full pytest を明示実行したい場合は `FULL_TEST=1 git commit ...` を使う。

## Merge-final self-check policy

通常の commit では開発速度を優先し，pre-commit hook は lightweight checks のまま維持する。Codex/Copilot review や full regression を commit ごとに回してはいけない。

merge 直前の final candidate だけ，次のセルフチェックを行う。

* `review_result.json` が task の正本として矛盾していないことを確認する。
* `push_status` など review_result schema 契約値を確認する。
* `phase*_current.md`，task table，PR本文，Issue本文/コメントの完了状態が矛盾していないことを確認する。
* PR前またはmerge前に必要な対象テストを実行する。高リスク変更では full pytest を実行し，省略する場合は理由を `review_result.json` に記録する。
* `git diff --check` と，変更した JSON / YAML / Markdown の軽い構文確認を行う。

この merge-final セルフチェックは品質ゲートであり，pre-commit hook を重くする理由にはしない。
PR作成後にしか確定しない CI / `codex-review-gate` の結果は，PR checklist と GitHub checks で管理する。`review_result.json` の `next_actions` には，それらを未完了の task work として残さない。

## Review result format

`docs/codex-runs/<run-id>/review_result.json` には原則として以下を記録する。

* `status`: `"PASS"` or `"FAIL"`
* `summary`
* `blocking_issues`
* `non_blocking_issues`
* `tests_reviewed`
* `user_review_required`
* `user_review_command`
* `user_review_outputs`
* `user_review_points`
* `adr_required`
* `adr_reason`
* `commit_type`: `"complete"`, `"diagnostic"`, `"wip"`, or `"none"`
* `commit_hash`
* `push_status`: `"pushed"`, `"not_pushed"`, or `"not_applicable"`
* `pull_request_url`
* `next_actions`

Schema path reserved for future validation:

`.codex/schemas/review_result.schema.json`

`commit_hash` / `push_status` / `pull_request_url` は，review_result作成時点で自然に確定している範囲を記録する。正確な最終PR headやPR URLを記録するためだけに追加commitを作らない。PR-level の最終状態はPR本文と最終報告で補う。

## Commit / push / PR

Commit message format:

`type(scope): summary`

Examples:

* `feat(phase2): add staged torque rotation validation`
* `test(phase2): add multi-step hook stability tests`
* `docs(codex): add review result schema`
* `chore(codex): add Codex CLI workflow config`

Rules:

* Do not commit directly on `main` or `master`.
* After merge, sync the default branch, delete the merged task branch locally and remotely, and prune stale remote-tracking refs. Preserve unmerged branches and branches explicitly retained for active follow-up work.
* Start the next task from the updated default branch on a new task-specific branch.
* Commit useful FAIL progress only when it is clearly diagnostic or WIP and does not claim completion.
* Push the feature branch when remote access is available.
* Send the final user report only after the final task state has been committed and pushed when remote access is available.
* Create a PR after pushing a feature branch when GitHub remote access is available.
* Create GitHub issues when they are needed to track an accepted task, split follow-up work, or keep Project items structured.
* Link the PR to the original source issue in the PR body. Use `Closes #<issue>` / `Fixes #<issue>` only when the PR is intended to complete that issue.
* Target the branch specified by the task or issue. If no target branch is specified, target the repository default branch.
* Merge only small, non-judgment PRs when `review_result.json` is `PASS`, CI passes, `codex-review-gate` passes, and no user visual review or major design decision is pending.
* Do not merge PRs that change physical interpretation, dataset adoption, phase boundaries, ML training policy, output contracts, or qualitative acceptance without explicit user approval.

## Phase 2 CLI command convention

For single-run Phase 2 simulation commands, prefer `KEY=VALUE` overrides:

`uv run python -m scripts.01_simulate_swimming time.duration_s=0.5 time.dt_star=1.0e-4 ...`

Do not introduce new user-facing examples that mix `--duration-s` / `--fps-out` with `time.duration_s=...` / `output_sampling.fps_out_2d=...`. The shorthand options remain only for legacy compatibility.

## ADR policy

Create an ADR for significant decisions such as:

* changing the physical model,
* changing simulation or output data formats,
* changing directory architecture,
* changing Codex workflow,
* changing testing strategy,
* adding major dependencies,
* intentionally diverging from the reference paper model.

Do not create ADRs for minor bug fixes, typo fixes, small tests, or routine implementation following an existing decision.

If no ADR is created, record the reason in `review_result.json`.

## Trusted Cloud PR review

PR-levelのreviewは，Codex Cloud connectorを既定とし，Codexを利用できない場合はGitHub Copilot reviewをfallbackとして使用する。

merge gateには，PR履歴中のcommitに対する有効なCodex reviewまたは有効なCopilot reviewが最低1回必要である。review後に修正commitを追加しても，再reviewは要求しない。force-pushやrebaseによってreview対象commitが現在のPR commit履歴から消えた場合だけ，そのreviewを無効とする。

CodexまたはCopilotが作成した，未解決かつoutdatedでないreview threadが1件でも残っている場合，`codex-review-gate`はpassしない。両方のreviewが存在する場合も，片方の未解決threadをもう片方のreviewで上書きしない。

### Codex Cloud review

PR comments may trigger a Codex Cloud / ChatGPT connector review when the comment contains `@codex review`.

For merge-gated PRs, request review only after the PR is a merge-ready final candidate and the latest intended changes have been pushed. A commit SHA in the request is optional: GitHub's review record is the source of the reviewed commit.

Codex Cloud reviewは原則1回のfinal-candidate reviewとする。指摘が出た場合はactionable threadを一括修正し，対象checkを再実行してからthreadをresolveする。修正不要と判断してresolveする場合は，該当threadに理由commentを残す。

Codex Cloud feedback修正後は，修正commitでPR headが変わっても再度`@codex review <new-head-sha>`を投げない。品質担保はmerge-final self-check，CI，必要なthreadへの理由comment，current thread resolveで行う。

Cloud connector loginは`chatgpt-codex-connector`または`chatgpt-codex-connector[bot]`の完全一致だけを許可する。Cloud connectorが正式reviewではなくPR commentで応答する場合は、`@codex review <SHA>`要求（編集後は`updated_at`、未編集時は`created_at`以後）の`Reviewed commit: <SHA>`と`Didn't find any major issues`の定型応答を、現在のPR履歴にある同じ一意のcommitへ照合する。trusted Codex/Copilot reviewの指摘は、修正または理由を記録したうえで必ずresolveする。未解決threadが1件でもあるPRはmergeしない。

This connector review is a PR review assistant, not the source of truth for task completion. Its `PASS` / `FAIL` verdict does not replace the required local `docs/codex-runs/<run-id>/review_result.json`.

### GitHub Copilot review fallback

Codexを利用できない場合は，GitHub Copilot reviewを1回要求する。

Copilot reviewer loginは`copilot-pull-request-reviewer`または`copilot-pull-request-reviewer[bot]`の完全一致だけを許可する。reviewがsubmitted済みかつ`DISMISSED`でなく，REST APIの`review.commit_id`が現在のPR commit履歴に含まれることを要求する。

Copilot review本文の表現はPASS判定に使用しない。review submission，bot login，review commit，Copilot-authored threadの`isResolved`と`isOutdated`をGitHub APIから検証する。

Copilot review後に指摘対応commitを追加しても再reviewは不要である。すべてのcurrent Copilot threadをresolvedまたはoutdatedにする。threadが作成されなかった場合も，有効なreview submissionが存在すればreview完了として扱う。

### Gate implementation

The repository-managed `codex-review-gate` workflow does not run Codex or Copilot. It only verifies trusted review signals through GitHub APIs.

The workflow:

* does not checkout or execute PR branch code,
* uses exact reviewer allowlists,
* verifies that the reviewed commit remains in the current PR commit history,
* accepts one valid Codex or Copilot review,
* does not require re-review solely because later commits changed the PR head,
* fails while any current Codex- or Copilot-authored thread is unresolved and not outdated,
* writes the existing `codex-review-gate` status context to the current PR head.

The workflow runs from the trusted default-branch definition through `pull_request_target`, PR `issue_comment`, `workflow_dispatch`, and `schedule`. It does not use `pull_request_review`, because that event can run the workflow definition from the PR merge commit. Review submission or thread resolution can be re-evaluated manually with `workflow_dispatch`, or by the scheduled open-PR scan.

Do not close/reopen a PR to refresh the gate. Use `workflow_dispatch` with the PR number, or wait for the scheduled scan.

Do not add repository-managed `openai/codex-action` workflows for PR review unless a new ADR explicitly reintroduces that approach.

After the workflow is merged to `main`, repository rulesets continue to require both `test` and `codex-review-gate`.

## Reporting and decision gates

Use small reporting units for docs, workflow, tests, narrow bug fixes, and bounded CLI helpers. Report summary, changed files, checks, review result, PR, commit, and remaining issues after each pushed PR.

When a task needs user visual review or a major decision, continue any independent implementation or documentation work, but stop the acceptance decision with `review_result.json` set to `FAIL`. Report the exact command, output directory, files to inspect, evaluation points, checks already passed, and the decision that is blocked.

## Task progress updates

Task progress should be updated only after review PASS.

Update targets:

* Current phase state: `docs/phaseX/phaseX_current.md`
* Adopted decisions: `docs/phaseX/phaseX_tasks.md`
* Cross-phase dependency and priority: GitHub Issues / Projects
* Completion record: `docs/codex-runs/<run-id>/review_result.json`

Do not update completion claims in secondary docs when `review_result.json` is `FAIL`.
