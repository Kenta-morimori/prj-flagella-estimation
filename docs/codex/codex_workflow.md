# Codex ワークフロー

この文書は、毎回読む必要のない Codex 運用の詳細をまとめる。

通常は `AGENTS.md` と対象 task の current doc だけを読み、完了条件・review_result・commit/push/PR・ADR・Cloud review の判断が必要なときだけこの文書を読む。

## 正本

`docs/codex-runs/<run-id>/review_result.json` は Codex task の完了状態に関する正本である。

task checkbox、commit message、PR 本文、Codex の final response は二次記録であり、`review_result.json` と矛盾してはいけない。

`review_result.json` の `PASS` は、ローカルの実装・文書・セルフチェックが完了していることを表す。PR 作成後の CI と trusted Cloud review（Codex または GitHub Copilot）は merge gate であり、PR checklist と GitHub checks で確認する。trusted review が未実施であることだけを理由に、ローカル完了済みの `review_result.json` を `FAIL` に戻さない。

## Issue の execution target

新規 Issue は GitHub Issue Form で Heavy/runtime execution target を必須選択する。コード実装と短時間 unit test は通常 Mac で行い、この target は長時間 simulation、sweep、render などの実行先を示す。

| Form value | Label | 許可する実行 |
| --- | --- | --- |
| `mac_only` | `execution:mac` | Mac local runtime のみ |
| `cs10_user_run` | `execution:cs10` | Mac 実装・短時間 check 後、User が cs10 で heavy job を実行 |
| `no_runtime` | `execution:none` | docs / review / workflow のみ |
| `triage_required` | `execution:triage` | read-only triage のみ |

独立 condition が8以上、または Mac 見積り wall time が30分超なら、`cs10_user_run` を選ぶ必須候補とする。Issue 作成・編集時の workflow が `execution:*` label を同期する。本文の target と label が不在・不一致、または `execution:triage` なら、Codex は実装・test・runtime を開始せず、target の triage を依頼する。既存 Issue は一括推測せず、着手時にこの Form 項目を追記して triage する。

`execution:cs10` かつ独立 condition が2以上なら、Codex は parallel job config、worker plan、`cs10_qualified`、dry-run、condition ごとの output 分離を確認するまで runtime を開始しない。serial 例外には、Issue runbook の具体的な技術理由と、開始前の User 明示承認 Issue コメント URL が必要である。

## Issue の Roadmap metadata

新規 Issue は Issue Form の必須 `Roadmap category (Milestone)` を選択する。`issue-roadmap-sync` workflow は Project #8 へ Issue を登録し、対応 Milestone と作成日（JST）の `Start date` を同期する。`Planned target date` は任意の `YYYY-MM-DD` 入力であり、指定時は `Target date` へ同期する。Issue close 時に Target date が未設定なら、その Issue の終了日（JST）を補完する。既に予定日があれば上書きしない。

この workflow は user-owned Project を更新するため、classic PAT（`repo` + `project` scope）を `PROJECT_AUTOMATION_TOKEN` として repository secret へ登録する。token が無い場合は明示的に失敗する。Issue Form 外で作成された Issue や不正な日付は `roadmap:triage`、reopen された Issue は `roadmap:needs-review` で明示し、metadata 修正まで実装に着手しない。

PR URL、最終 PR head SHA、push 後の状態など、PR 作成後にしか確定しない動的情報を tracked `review_result.json` へ後追い同期するためだけの commit は作らない。これらは PR 本文、GitHub checks、最終ユーザー報告に記録する。

file-changing な Issue 実装では、初回完了報告は PR 作成後まで送らない。local PASS、commit、push、source Issue を参照する PR 作成の後に、PR URL と未完了の merge gate を報告する。commentary の進捗共有は可能だが、PR 作成前に実装完了・成果物・PR 候補として報告してはならない。例外はユーザーが明示的に PR 不要とした場合、または PR 作成が失敗した場合だけであり、後者は試行内容と concrete blocker を報告する。

## Run ID

形式:

`YYYYMMDD_HHMMSS_<phase>_<task-id>`

例:

`docs/codex-runs/20260530_142233_phase2_0037/review_result.json`

## 完了ポリシー

`PASS` の完了には以下が必要である。

1. 要求された実装または文書変更が完了している。
2. relevant tests/checks が `PASS`、または未実行理由が明確である。
3. local review step が完了している。PR 作成後の CI / trusted Cloud review は merge gate として別管理する。
4. `docs/codex-runs/<run-id>/review_result.json` が `"status": "PASS"` である。
5. work log / review result が保存されている。
6. final state が commit 済みである。
7. remote access があれば push 済みである。
8. pushed feature branch なら PR が作成済みである。

`FAIL` result は完了ではない。ただし Phase 2 では、有用な診断進捗を `diagnostic`、`wip`、`docs`、`test` 相当の commit として保存してよい。

有用な `FAIL` の例:

* collapse / fly-away / hook drift / no_bundle 条件を再現した。
* failing test で次の target behavior を定義した。
* 物理モデル差分や数値上の不一致を記録した。
* 部分実装で原因範囲を狭めた。

## Test ポリシー

pre-commit hook の既定は lightweight checks とする。

既定コマンド:

* `uv run ruff format --check .`
* `uv run ruff check .`
* `uv run pytest -q -m light`

`light` は commit 時に固定実行しても負担が小さい、短時間・deterministic・library-level の test を指す。初期運用では対象を広げすぎず、明らかに軽い test だけを明示的に marker 付与する。

full pytest は削除しない。以下では `uv run pytest -q` を実行する。

* 物理モデル、geometry、hook、flagella、body、torque、force、potential、hydrodynamics の変更時
* simulation core の変更時
* output format、manifest、CSV schema の変更時
* dataset 生成仕様の変更時
* PR 作成前または merge 前
* GitHub Actions CI

docs-only、planning-only、workflow-only 変更では full pytest を既定要求しない。ただし未実行の場合も、必要なら `review_result.json` に理由を書く。

hook で full pytest を明示実行したい場合は `FULL_TEST=1 git commit ...` を使う。

## Merge 前の最終セルフチェックポリシー

通常の commit では開発速度を優先し、pre-commit hook は lightweight checks のまま維持する。Codex/Copilot review や full regression を commit ごとに実行してはいけない。

merge 直前の final candidate だけ、次のセルフチェックを行う。

* `review_result.json` が task の正本として矛盾していないことを確認する。
* `push_status` など review_result schema 契約値を確認する。
* `phase*_current.md`、task table、PR 本文、Issue 本文/コメントの完了状態が矛盾していないことを確認する。
* PR 前または merge 前に必要な対象 test を実行する。高リスク変更では full pytest を実行し、省略する場合は理由を `review_result.json` に記録する。
* `git diff --check` と、変更した JSON / YAML / Markdown の軽い構文確認を行う。

この merge 前の最終セルフチェックは品質ゲートであり、pre-commit hook を重くする理由にはしない。PR 作成後にしか確定しない CI / `codex-review-gate` の結果は、PR checklist と GitHub checks で管理する。`review_result.json` の `next_actions` には、それらを未完了の task work として残さない。

## Review result の形式

`docs/codex-runs/<run-id>/review_result.json` には原則として以下を記録する。

* `status`: `"PASS"` または `"FAIL"`
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
* `commit_type`: `"complete"`、`"diagnostic"`、`"wip"`、または `"none"`
* `commit_hash`
* `push_status`: `"pushed"`、`"not_pushed"`、または `"not_applicable"`
* `pull_request_url`
* `next_actions`

将来の validation 用 schema path:

`.codex/schemas/review_result.schema.json`

`commit_hash` / `push_status` / `pull_request_url` は、review_result 作成時点で自然に確定している範囲を記録する。正確な最終 PR head や PR URL を記録するためだけに追加 commit を作らない。PR-level の最終状態は PR 本文と最終報告で補う。

## Commit / push / PR

commit message の形式:

`type(scope): summary`

例:

* `feat(phase2): add staged torque rotation validation`
* `test(phase2): add multi-step hook stability tests`
* `docs(codex): add review result schema`
* `chore(codex): add Codex CLI workflow config`

ルール:

* `main` または `master` に直接 commit しない。
* merge 後は default branch を同期し、merge 済み task branch を local と remote で削除し、古い remote-tracking ref を prune する。merge 前の branch と、継続作業のため明示的に保持する branch は保全する。
* 次の task は、更新済み default branch から新しい task-specific branch で開始する。
* 有用な `FAIL` 進捗は、明確に diagnostic または WIP であり完了を主張しない場合だけ commit する。
* remote access があれば feature branch を push する。
* GitHub remote access があれば、feature branch を push した後に PR を作成する。
* 初回完了報告は、final task state の commit・push と source-Issue-linked PR の作成後にだけ送る。例外は、PR を作成しないという明示的なユーザー指示、または PR 作成失敗の試行と concrete blocker の報告だけである。
* 新規 GitHub Issue / sub-issue は、受入済み task の追跡、follow-up の分割、Project 構造の維持に必要な場合だけ作成する。残作業の候補を見つけただけでは作成せず、具体案を提示してユーザー承認を得る。
* PR 本文で元の source Issue に PR を紐付ける。PR がその Issue を完了する意図のときだけ `Closes #<issue>` / `Fixes #<issue>` を使う。
* task または Issue が指定した branch を target とする。target branch の指定がなければ repository default branch を target とする。
* `review_result.json` が `PASS`、CI が pass、`codex-review-gate` が pass であり、ユーザー visual review や major design decision が未解決でない場合だけ、小規模で判断不要な PR を merge する。
* 物理解釈、dataset adoption、phase boundary、ML training policy、output contract、または qualitative acceptance を変更する PR は、明示的なユーザー承認なしに merge しない。

## PRマージ後のIssue引継ぎ

source Issue に紐づく PR が merge されたら、PR、source Issue、ローカル `review_result.json` を照合してから完了状態を報告する。次を確認する。

1. merge された PR と source Issue の対応、`Closes` / `Fixes` の有無、source Issue の open / closed 状態。
2. Issue の全受入条件、PR の実装範囲、relevant check、`review_result.json: PASS` が一致すること。
3. すべて満たす場合だけ、必要な Issue checkbox、状態、完了記録を更新する。未達の条件があれば、完了として扱わない。
4. 継続親Issueは、子 PR または child Issue が merge / close しただけでは閉じない。親自身の受入条件と継続目的を別に確認する。
5. 残作業があれば、範囲、受入条件、`Heavy/runtime execution target` を含む具体的な後続 task を source Issue の更新と最終報告で提示する。新規 Issue / sub-issue の作成は、ユーザー承認後に行う。

GitHub 更新の権限が task にない場合は、必要な Issue 更新内容と後続 task 案を最終報告に明記し、外部状態を推測して更新済みとしない。

## Phase 2 CLI command の慣例

単一 run の Phase 2 simulation command では、`KEY=VALUE` override を優先する。

`uv run python -m scripts.01_simulate_swimming time.duration_s=0.5 time.dt_star=1.0e-4 ...`

`--duration-s` / `--fps-out` と `time.duration_s=...` / `output_sampling.fps_out_2d=...` を混在させる新しい user-facing example を導入しない。shorthand option は legacy compatibility のためだけに残す。

## ADR ポリシー

次のような重要な判断には ADR を作成する。

* 物理モデルの変更
* simulation または output data format の変更
* directory architecture の変更
* Codex workflow の変更
* testing strategy の変更
* major dependency の追加
* 参照論文モデルから意図的に乖離する変更

軽微な bug fix、typo fix、小規模 test、既存判断に従う通常実装には ADR を作成しない。

ADR を作成しない場合は、理由を `review_result.json` に記録する。

## 信頼できる Cloud PR review

PR-level の review は、Codex Cloud connector を既定とし、Codex を利用できない場合は GitHub Copilot review を fallback として使用する。

merge gate には、PR 履歴中の commit に対する有効な Codex review または有効な Copilot review が最低1回必要である。review 後に修正 commit を追加しても、re-review は要求しない。force-push や rebase によって review 対象 commit が現在の PR commit 履歴から消えた場合だけ、その review を無効とする。

Codex または Copilot が作成した、未解決かつ outdated でない review thread が1件でも残っている場合、`codex-review-gate` は pass しない。両方の review が存在する場合も、片方の未解決 thread をもう片方の review で上書きしない。

### Codex Cloud review

PR comment に `@codex review` を含めると、Codex Cloud / ChatGPT connector review を trigger できる。

merge-gated PR では、PR が merge-ready final candidate となり、意図した最新変更が push された後にだけ review を依頼する。依頼に commit SHA を含めるかは任意である。GitHub の review record を review 対象 commit の正本とする。

Codex Cloud review は原則1回の final-candidate review とする。指摘が出た場合は actionable thread を一括修正し、対象 check を再実行してから thread を resolve する。修正不要と判断して resolve する場合は、該当 thread に理由 comment を残す。

Codex Cloud feedback 修正後は、修正 commit で PR head が変わっても再度 `@codex review <new-head-sha>` を投稿しない。品質担保は merge 前の最終セルフチェック、CI、必要な thread への理由 comment、current thread の resolve で行う。

Cloud connector login は `chatgpt-codex-connector` または `chatgpt-codex-connector[bot]` の完全一致だけを許可する。Cloud connector が正式 review ではなく PR comment で応答する場合は、`@codex review <SHA>` 要求（編集後は `updated_at`、未編集時は `created_at` 以後）の `Reviewed commit: <SHA>` と `Didn't find any major issues` の定型応答を、現在の PR 履歴にある同じ一意の commit へ照合する。trusted Codex/Copilot review の指摘は、修正または理由を記録したうえで必ず resolve する。未解決 thread が1件でもある PR は merge しない。

この connector review は PR review assistant であり、task completion の正本ではない。その `PASS` / `FAIL` verdict は、必要なローカル `docs/codex-runs/<run-id>/review_result.json` を置き換えない。

### GitHub Copilot review fallback

Codex を利用できない場合は、GitHub Copilot review を1回要求する。

Copilot reviewer login は `copilot-pull-request-reviewer` または `copilot-pull-request-reviewer[bot]` の完全一致だけを許可する。review が submitted 済みかつ `DISMISSED` でなく、REST API の `review.commit_id` が現在の PR commit 履歴に含まれることを要求する。

Copilot review 本文の表現は `PASS` 判定に使用しない。review submission、bot login、review commit、Copilot-authored thread の `isResolved` と `isOutdated` を GitHub API から検証する。

Copilot review 後に指摘対応 commit を追加しても re-review は不要である。すべての current Copilot thread を resolved または outdated にする。thread が作成されなかった場合も、有効な review submission が存在すれば review 完了として扱う。

### Gate 実装

repository-managed `codex-review-gate` workflow は Codex や Copilot を実行しない。GitHub API を通じて trusted review signal だけを検証する。

この workflow は次を満たす。

* PR branch code を checkout または実行しない。
* reviewer allowlist を完全一致で使用する。
* review 対象 commit が現在の PR commit 履歴に残ることを検証する。
* 有効な Codex または Copilot review を1件受け入れる。
* 後続 commit により PR head が変わったことだけを理由に re-review を要求しない。
* current Codex- または Copilot-authored thread が unresolved かつ outdated でない間は fail する。
* 既存の `codex-review-gate` status context を現在の PR head に書き込む。

この workflow は、trusted default-branch definition から `pull_request_target`、PR `issue_comment`、`workflow_dispatch`、`schedule` で実行する。`pull_request_review` event は PR merge commit の workflow definition を実行し得るため使用しない。review submission または thread resolution は、`workflow_dispatch` か scheduled open-PR scan により手動で再評価できる。

gate を更新するために PR を close/reopen しない。PR number を指定した `workflow_dispatch` を使うか、scheduled scan を待つ。

新しい ADR でこの方針を明示的に再導入しない限り、PR review のための repository-managed `openai/codex-action` workflow を追加しない。

workflow が `main` に merge された後も、repository ruleset は `test` と `codex-review-gate` の両方を要求し続ける。

## 報告と判断ゲート

docs、workflow、test、限定的な bug fix、範囲を限定した CLI helper には小さな報告単位を使う。file-changing な Issue 作業では、source-Issue-linked PR が存在した後にだけ、summary、changed files、checks、review result、PR、commit、remaining issues を報告する。進捗 commentary はその前でも許可される。PR 作成失敗は concrete blocker としてのみ報告してよい。

task がユーザー visual review または major decision を必要とする場合、独立した実装または文書作業は継続してよいが、受入判断は `review_result.json` を `FAIL` にして止める。exact command、output directory、確認対象 file、evaluation point、すでに pass した check、block されている decision を報告する。

## Task の進捗更新

task の進捗は review `PASS` 後にだけ更新する。

更新対象:

* 現在の phase state: `docs/phaseX/phaseX_current.md`
* 採用済み decision: `docs/phaseX/phaseX_tasks.md`
* cross-phase dependency と priority: GitHub Issues / Projects
* completion record: `docs/codex-runs/<run-id>/review_result.json`

`review_result.json` が `FAIL` のとき、二次文書の完了主張を更新しない。
