# Context Routing

このreferenceは，Issue単位作業で読む文書を絞るためのrouting表である．

## Always Check First

- `git status --short --branch`
- user requestの最新内容
- 対象Issue / PRの本文と最新コメント
- 対象Phase

## Repository Policy

読む条件:

- repository全体の禁止事項を確認する
- branch，commit，push，mergeの扱いを判断する
- testing，review，documentationの共通規則を確認する

読む文書:

- `AGENTS.md`

注意:

- `AGENTS.md`はrepository policyの正本である．
- skillや他docsと矛盾した場合は`AGENTS.md`を優先する．
- 詳細なPhase文書規則は`docs/codex/phase_document_policy.md`を正本とする．

## Current Phase State

読む条件:

- Phaseの現在地を確認する
- current baselineを確認する
- Active work，Next queue，blockerを確認する
- 今回読むべき文書を確認する

読む文書:

- `docs/phaseX/phaseX_current.md`

読み方:

- 最初に対象Phaseのcurrentを読む．
- currentから必要なcontract，guide，validationへ進む．
- project-level docsは全体地図やPhase間依存が必要な場合だけ読む．
- 完了履歴を探すためにcurrentを使用しない．

## Phase-Specific Guide

読む条件:

- Phase固有の恒常的な実行規則が必要
- CLI，input / output boundary，長時間実行方針を確認する
- モデル，dataset，trainingの解釈上の注意を確認する
- Phase固有の目視レビュー観点を確認する

読む文書:

- `docs/phaseX/phaseX_guide.md`

注意:

- guideが存在しないPhaseでは省略する．
- guideからActive Issue，個別run結果，temporary blockerを探さない．

## Past Decisions

読む条件:

- 過去に何を採用，不採用，置換，保留としたか確認する
- モデルや実行条件を変更した背景を確認する
- dataset，schema，pipeline，training policyの採択理由を確認する
- 過去結果の解釈を再利用する
- 既存判断を変更してよいか確認する

読む文書:

- `docs/phaseX/phaseX_tasks.md`

探し方:

```bash
rg -n \
  "Issue #<num>|P<X>-|<decision-id>|<model-name>|<dataset-version>|<keyword>" \
  docs/phaseX/phaseX_tasks.md
```

tasksから確認する情報:

- status
- background
- change
- comparison
- result
- interpretation
- decision
- evidence

tasksから原則として確認しない情報:

- branch履歴
- 完全なacceptance criteria
- 全実行command
- output directory一覧
- 実装ファイル一覧
- 時系列の作業日誌

これらはIssue，PR，Git履歴を参照する．

## Live Contracts And Implementation Sources

読む条件:

- 現在の実装契約を確認する
- schemaやfield名を確認する
- config defaultやparameter値を確認する
- validation gateを確認する
- feature定義やPhase handoffを確認する
- active validationを実行する

読む対象:

- currentまたはguideから参照されるtask-specific doc
- `schemas/`
- `conf/`
- code
- tests
- active validation docs
- dataset registry
- handoff contract

正本の優先:

- parameter値: config
- machine-readable field: schema
- gate実装: codeとtest
- detailed design reason: ADR
- 採択判断: tasks

docsとmachine-readable sourceが矛盾する場合は，矛盾を報告し，無断で一方へ合わせない．

## Phase Document Maintenance

読む条件:

- currentまたはtasksを再構成する
- task-specific docを統合・削除する
- Phase文書間の重複を解消する
- Phase 2の構成を他Phaseへ適用する
- 文書削除前の情報移行を確認する

読む順序:

1. `docs/codex/phase_document_policy.md`
2. `tools/codex/skills/phase-document-maintenance/SKILL.md`
3. `tools/codex/skills/phase-document-maintenance/references/consolidation-checklist.md`
4. 対象Phaseのcurrentとtasks
5. 統合候補文書
6. 必要なADR，Issue，PR，config，schema，test

文書整理時も，対象Phaseの全ファイルを無条件に全文読込しない．
まずファイル一覧と`rg -n`で候補を絞る．

## Issue And PR Context

読む条件:

- scopeとacceptance criteriaを確認する
- branchやPRの完了範囲を確認する
- tasksに記録されていない詳細結果を確認する
- 過去判断のevidenceを確認する

読む対象:

- source Issue
- target PR
- related Issue / PR
- review comments

Issue / PRの情報をtasksへそのまま複製せず，判断に必要な内容だけを圧縮する．

## ADR

読む条件:

- 重要な設計理由を確認する
- reference modelとの差分を確認する
- architecture，schema，dataset policy，Phase boundaryの変更を確認する

読み方:

- tasksから参照されたADRを優先する．
- 結論だけで十分な場合はtasksの要約を使用する．
- ADR本文をcurrentやtasksへ再複製しない．

## Codex Run Logs

読む条件:

- 過去runのPASS / FAILを確認する
- commit hash，PR，blocking issue，next actionを確認する
- 既存branchの作業を引き継ぐ

探し方:

```bash
rg -n \
  "source issue|pull_request_url|commit_hash|blocking_issues|next_actions|Issue #<num>" \
  docs/codex-runs
```

読み方:

1. `review_result.json`
2. 必要な場合だけ同じrunの`work_log.md`
3. summaryで不足する場合だけraw logsやgenerated output

長いlogs，CSV，generated outputはデフォルトで全文表示しない．

## Removed Or Historical Context

現行repositoryでは`prompts/`をsource of truthとして扱わない．

削除済み文書や古いacceptance criteriaの詳細が必要な場合は，以下の順で確認する．

1. tasks
2. ADR
3. Issue / PR
4. `review_result.json`
5. Git history
6. `work_log.md`
7. generated output

削除済み文書を通常routingへ再追加しない．
