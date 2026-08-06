---
name: phase-document-maintenance
description: Use when auditing, restructuring, consolidating, migrating, or deleting docs/phase*/ documents, including phase current files, tasks decision records, guides, contracts, validation documents, plans, and reports.
---

# Phase Document Maintenance

## Overview

このskillは，`docs/phase*/`の文書を，後続開発に必要な判断を失わずに整理・軽量化するためのworkflowである．

文書構成と保持基準の正本は`docs/codex/phase_document_policy.md`とする．このskillでは，policyを実際の調査・統合・削除作業へ適用する手順を定める．

対象はPhase 2に限定せず，Phase 3，Phase 4以降にも適用する．

## When To Use

以下の場合に使用する．

- `phaseX_current.md`を現在地だけに整理する
- `phaseX_tasks.md`を判断記録として再構成する
- 完了済みplan，design note，reportをtasksへ統合する
- task-specific docの保持・削除を判断する
- Phase文書間の重複や矛盾を解消する
- Phase 2で確立した構成をPhase 3/4へ展開する
- semantic decisionを伴うタスク完了時に文書を更新する
- 文書削除前に情報移行と参照整合性を確認する

単純なtypo，表記修正，リンク修正だけの場合は，このskillを全面的に適用する必要はない．

## Required Context

作業開始時は，必要な範囲に限定して以下を読む．

1. `AGENTS.md`
2. user requestと対象Issue / PR
3. `docs/codex/phase_document_policy.md`
4. 対象Phaseの`phaseX_current.md`
5. 対象Phaseの`phaseX_tasks.md`の関連section
6. 統合・削除候補のtask-specific docs
7. 必要なADR，config，schema，test，Issue / PR

長い文書は最初から全文を読まず，`rg -n`でIssue番号，task ID，model名，dataset version，config keyを検索する．

## Scope First

最初に以下を明確にする．

- 対象Phase
- 対象Issue / PR
- auditのみか，実際に変更するか
- 文書削除を許可されているか
- semanticな仕様変更を含むか
- current，tasks，guideの更新範囲
- Phase 3/4などのscope外

文書整理のみがscopeの場合，既存のmodel，parameter，dataset採択，schema，training policyを変更しない．

## Workflow

### 1. Inventory

対象Phaseの文書を列挙し，以下へ分類する．

- `current`
- `tasks`
- `guide`
- `live contract`
- `active validation`
- `reusable report`
- `obsolete history`
- `unclear`

ファイル名だけで分類しない．本文，参照元，config，schema，test，ADR，Issue / PRを確認する．

詳細な確認項目とinventory形式は`references/consolidation-checklist.md`を使用する．

### 2. Extract Decisions

統合候補から，後続開発で必要な情報を抽出する．

- background
- change
- comparison conditions
- decision-bearing result
- interpretation and limitations
- adopted / rejected / replaced / pending decision
- evidence

単なる実行履歴と判断情報を分ける．

### 3. Consolidate Tasks

`phaseX_tasks.md`はIssue単位のwork logではなく，判断単位で整理する．

複数Issue，PR，runが同じ判断へ収束した場合は，1つのdecision sectionへ統合してよい．

entry形式，status，数値の保持基準は`docs/codex/phase_document_policy.md`に従う．

以下は通常tasksへ移さない．

- branch名
- 完全なacceptance criteria
- 全実行コマンド
- output directory一覧
- 実行日時
- 全seed / sweep cell
- PR本文と重複する実装一覧

### 4. Update Current

`phaseX_current.md`には，現在の作業開始に必要な情報だけを残す．

- current baseline
- Active work
- Next queue
- blockers
- context routing

完了済みタスク，長い結果表，過去のoutput pathを残さない．

### 5. Preserve Live Documents

以下は原則として個別文書を維持する．

- schema / output-column definition
- feature definition
- validation contract
- parameter mapping
- dataset version registry
- Phase handoff contract
- active validation procedure
- 後続解析で再利用する大規模report

tasksへ圧縮すると実装情報が欠損する場合は，文書を削除しない．

### 6. Delete Obsolete History

完了済みplan，temporary note，候補比較，実行報告は，以下を満たす場合に削除できる．

- decision-bearing informationがtasksまたはADRへ移行済み
- live contract情報が別の正本に残る
- 詳細なevidenceを辿れる
- 参照元を更新できる
- 採用状態がcurrent / tasks / ADRと一致する

不明な場合は`unclear`または`active validation`として保持し，削除を保留する．

### 7. Update References

削除・移動・名称変更後は，以下を対象に旧path，旧filename，deprecated terminologyを検索する．

- `AGENTS.md`
- `README.md`
- `docs/`
- `tools/`
- `scripts/`
- `conf/`
- `schemas/`
- `tests/`

参照先は，原則としてtasks，ADR，live contract，Issue / PRのいずれかへ置き換える．

### 8. Verify Consistency

最低限，以下を確認する．

- currentとtasksの採用状態が一致する
- tasksとADRの結論が一致する
- parameter値はconfigと一致する
- field名はschemaと一致する
- gate値はcode / testと一致する
- 削除済み文書への参照が残っていない
- currentが履歴で再肥大化していない
- tasksが実行履歴の複製になっていない

### 9. Complete The Issue Workflow

文書整理後のreview，`review_result.json`，commit，push，PRは`flagella-issue-workflow`とその`references/completion-policy.md`に従う．

このskillは，Issue作業全体の完了条件を置き換えない．

## Safety Rules

- 情報移行前に文書を削除しない．
- ファイル名だけで削除判断しない．
- active validationの採否を推測しない．
- user decisionがpendingの結果を`adopted`にしない．
- Git履歴に残ることだけを理由に，通常開発で必要な判断を省略しない．
- 文書整理だけのtaskでsemanticなmodelやdataset条件を変更しない．
- evidenceが不足する場合は断定せず，保留理由を記録する．

## Verification

最低限:

```bash
git diff --check
```

削除・名称変更時:

```bash
rg -n "<deleted-path>|<deleted-file-name>" \
  AGENTS.md README.md docs tools scripts conf schemas tests
```

skill追加・更新時:

- YAML frontmatterが有効
- `name`とdirectory名が一致
- 記載したpathが存在
- policyおよび他skillと矛盾しない

文書のみの変更ではfull pytestを必須としない．code，config，schema，testを変更した場合は，対応するtargeted testを実行する．

## Completion Report

最終報告には以下を含める．

- 更新したcurrent
- 更新したtasks
- 統合した判断
- 維持したlive documents
- 削除したobsolete documents
- 削除を保留したdocumentsと理由
- 更新した参照
- 実行したcheck
- semantic changeの有無
- review result，commit，push，PR状態
