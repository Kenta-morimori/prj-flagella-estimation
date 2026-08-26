---
name: flagella-issue-workflow
description: Use for issue-driven Codex work in this repository, including planning, implementation, diagnostics, testing, result analysis, review_result logging, commits, pushes, PR creation, and context-efficient routing of project documents.
---

# Flagella Issue Workflow

## Overview

このskillは，`prj-flagella-estimation`のIssue単位作業を，必要な文書だけを読んで進めるためのworkflowである．

目的は，`AGENTS.md`の規則を守りながら，方針検討，実装，実行，分析，文書更新，完了処理を一貫して行い，contextと実行コストを抑えることである．

Phase文書の構成変更，統合，移行，削除を行う場合は，`phase-document-maintenance` skillを併用する．

## Start

1. `git status --short --branch`で現在branchと作業treeを確認する．
2. user request，対象Issue，対象PRを特定する．
3. Issue Formの`Heavy/runtime execution target`、condition数、Mac wall time見積り、`execution:*` labelを確認する。不在・不一致・`execution:triage`ならread-only triage以外を開始しない．
4. 対象Phaseを特定する．
5. `main`または`master`上で直接作業しない．
6. 作業タイプを分類する．
7. semantic decisionを伴うか確認する．

作業タイプ:

- `planning`: 方針検討，task decomposition，Issue draft，実装前整理
- `implementation`: code，config，schema，test，docsの変更
- `diagnostic`: simulation，sweep，失敗条件保存，原因切り分け
- `review-only`: 差分レビュー，CI確認，既存結果の整理
- `workflow`: Codex運用，skill，review result，PR方針の変更
- `documentation-maintenance`: Phase文書の監査，統合，移行，削除

## Context Routing

常に全資料を読まず，作業タイプと対象Phaseに応じて読む文書を絞る．

詳細なroutingは`references/context-routing.md`を参照する．

基本順序:

1. `AGENTS.md`
2. user requestと対象Issue / PR
3. `docs/phaseX/phaseX_current.md`
4. `docs/phaseX/phaseX_guide.md`が存在し，Phase固有規則が必要な場合
5. `docs/phaseX/phaseX_tasks.md`の関連section
6. live schema，contract，config，test，active validation
7. ADR
8. Issue / PR履歴，Git履歴，`review_result.json`

追加routing:

- 完了条件，`review_result.json`，commit，push，PRについては`references/completion-policy.md`を参照する．
- contextと実行量の削減については`references/resource-reduction.md`を参照する．
- Phase文書の整理・統合・削除では，`phase-document-maintenance` skillと`docs/codex/phase_document_policy.md`を参照する．
- Codex workflow全体の詳細が必要な場合のみ`docs/codex/codex_workflow.md`を参照する．

## Workflow

### 1. 方針検討

- Issue本文，current，必要なtasks section，live contractだけでscopeを決める．
- acceptance criteria，必要なtest / simulation，目視レビュー要否を分ける．
- 過去判断が必要な場合は，tasksをIssue番号，decision ID，model名，dataset version，関連keywordで検索する．
- tasksの要約だけで理由が不足する場合に限り，ADR，Issue，PR，Git履歴へ進む．
- planningのみの依頼では，重いsimulationやfull pytestを実行しない．

### 2. 実装

- 実装開始時にexecution target、condition数、Mac見積り、許可される実行をユーザーへ短く報告する．
- `execution:mac`はMac local、`execution:none`はruntimeなし、`execution:cs10`はMac実装・短時間check後にUser-run cs10 parallel jobへ進む。`execution:triage`はread-only investigationに限定する．
- 変更を依頼scope内に限定する．
- `scripts/`はuser-facing orchestration，`src/`は再利用可能な実装，`tools/codex/`はCodex workflow補助に使用する．
- config，schema，testなど既存の正本を確認し，同じ情報をdocsへ複製しない．
- Phase 2の物理モデル変更では，reference model，repository implementation，numerical stabilization，project-specific extensionを区別する．
- semantic decisionが発生した場合は，完了時にcurrent，tasks，ADRの更新要否を確認する．

### 3. 実行

- 最小のtargeted testから開始する．
- 重い処理は，short representative，targeted test，sweep，full executionの順で段階的に行う．
- 長時間simulation，sweep，training，renderは，ユーザーから明示的に依頼されない限り実行しない．
- user executionとする場合は，command，expected output，evaluation points，実行済みcheckを提示する．

### 4. 結果分析

- PASS / FAILだけでなく，blocking issue，non-blocking issue，次の切り分け対象を記録する．
- 意思決定に影響した結果と，補助的なdiagnostic resultを分ける．
- 結果が既存判断を支持するか，置換するか，保留するかを明示する．
- Phase 2のcollapse，fly-away，hook drift，no bundleなどは，再現条件が有用であればdiagnostic progressとして扱う．

### 5. Phase文書更新

semantic decisionを伴う場合は，以下を確認する．

- `phaseX_current.md`の更新
- `phaseX_tasks.md`の更新
- ADRの作成または更新
- task-specific docの保持・削除
- stale referenceの更新

Phase文書の再構成，情報移行，削除を伴う場合は，このskill内で詳細手順を重複させず，`phase-document-maintenance` skillへ委譲する．

### 6. 報告・完了

- `references/completion-policy.md`に従って`review_result.json`を作成する．
- PASS完了時はcommit，push，PRまで行う．
- FAILでも有用な診断結果は，diagnostic，wip，docs，test相当のcommitとして保存できる．
- FAILを完了扱いにしない．
- 文書変更時は，current，tasks，ADR，維持文書，削除文書を最終報告に含める．

## Resource Discipline

- `SKILL.md`にproject historyやPhase固有の詳細を蓄積しない．
- 長い規則やchecklistはreferenceまたはpolicyへ分離する．
- 既存docsは`rg -n`で候補を絞ってから読む．
- 長いMarkdown，logs，CSV，generated outputをデフォルトで全文表示しない．
- 過去runは`review_result.json`を先に読み，必要な場合だけ`work_log.md`を読む．
- currentへ完了履歴を追記し続けない．
- tasksへIssue，PR，command，output pathの全文をコピーしない．
- full test，long simulation，video renderは完了条件に必要な場合だけ実行する．
- 最終報告では，変更内容，実行したcheck，未実行項目，残Issueを簡潔に示す．
