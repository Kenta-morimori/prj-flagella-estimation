# Resource Reduction

このreferenceは，Codexのcontext使用量，tool実行量，simulation / test実行コストを抑えるための方針を定める．

文書構成と保持基準は`docs/codex/phase_document_policy.md`を正本とする．

## Principles

削減対象:

- 毎回読み込む文書量
- 同じ情報の重複
- 不要なIssue / PR / Git履歴の探索
- 長いMarkdown，CSV，log，generated outputの全文読込
- 不要なfull test
- 不要なlong simulation，sweep，training，render
- 完了済みplanやtemporary reportの蓄積

削減してはいけないもの:

- 判断の背景
- model，condition，schema，pipelineの変更内容
- 意思決定に影響した結果
- 結果の解釈
- 採用，不採用，置換，保留
- 再確認可能なevidence
- 現在利用するschema，contract，config，test
- reproducibilityに必要なraw output

## Progressive Context

基本的な読込順序:

1. `AGENTS.md`
2. user request
3. target Issue / PR
4. target Phaseのcurrent
5. Phase guide
6. tasksの関連section
7. live contract，config，schema，test
8. ADR
9. Issue / PR history
10. review result
11. Git history，work log，raw output

前段で十分な場合は，後段を読まない．

`phaseX_tasks.md`は全文を毎回読まず，Issue番号，decision ID，model名，dataset version，keywordで検索する．

## Phase Documents

### Current

`phaseX_current.md`には現在地だけを保持する．

- 完了履歴を追加し続けない
- 全open Issueを列挙しない
- 全metricや全output pathを記載しない
- config値を複製しない
- 詳細結果はtasks，ADR，reportへ分離する

### Tasks

`phaseX_tasks.md`は判断単位で圧縮する．

残すもの:

- background
- change
- representative comparison
- decision-bearing result
- interpretation
- decision
- evidence

残さないもの:

- Issue / PR本文の複製
- 全acceptance criteria
- 全command
- 全seed，sweep，time-step結果
- branch履歴
- output path一覧
- 時系列work log

### Task-Specific Documents

新規文書を作成する前に確認する．

- tasks entryだけで十分ではないか
- 既存guideまたはcontractへ追記できないか
- ADRとして残すべきか
- Issue / PRだけで十分ではないか

完了済みplan，temporary memo，comparison reportは，判断情報をtasksまたはADRへ移行した後に保持要否を再評価する．

文書の統合・削除は`phase-document-maintenance` skillを使用する．

## Lightweight And Heavyweight Work

軽量:

- Issue整理
- planning
- task decomposition
- docs proposal
- current / tasks更新
- Phase文書の監査・統合
- skillやworkflow更新
- targeted testだけで十分な小規模変更

重量:

- simulation sweep
- long simulation
- 3D / 2D render
- dataset生成
- ML training
- full pytest
- physics，schema，共通output contractの変更

運用:

1. static check
2. targeted test
3. short representative run
4. limited sweep
5. full testまたはlong execution

の順で段階的に進める．

docs-only，planning-only，workflow-only変更ではfull pytestを必須としない．

## Large Files

長いMarkdown，logs，CSV，generated outputは，必要な範囲だけ読む．

優先順位:

1. current
2. tasks section
3. compact summary
4. manifest
5. review result
6. selected columns / rows
7. raw file

避ける操作:

- 長いMarkdownの無条件な全文表示
- wide CSVへの無制限な`cat`
- 全列を含む`head` / `tail`
- 全run logの横断的な全文読込
- generated outputを最初の情報源として扱う

raw outputはreproducibilityのため保持してよいが，Codexの日常分析ではcompact summaryを優先する．

## Issue Scope

Issueまたはuser requestに以下があると，contextと実行量を抑えやすい．

```markdown
## Task type

- planning / implementation / diagnostic / review-only / workflow / documentation-maintenance

## Scope

- Target phase:
- Target files/modules:
- Out of scope:

## Execution budget

- Tests:
- Simulation:
- Render:
- Training:
- Full pytest required:

## Documentation impact

- current update:
- tasks update:
- ADR:
- document deletion:

## Completion target

- PASS completion required:
- Diagnostic FAIL can be committed:
- User review required:
```

すべての項目を毎回必須にはしない．
対象Phase，変更範囲，重い処理の実行可否が分かればよい．

## Planning And Execution

方針が不確かな重量taskでは，planningとimplementationを分離する．

分離が有効な例:

- 物理モデル候補が複数ある
- dataset採択条件が未確定
- schema変更範囲が不明
- long simulation前に短時間評価が必要
- user visual reviewが必要

小規模なdocs整理やworkflow修正は，同一PRでplanningから実装まで行ってよい．

## Completion Efficiency

- 過去runは`review_result.json`を先に確認する．
- `work_log.md`は要約で不足する場合だけ読む．
- semantic decisionを伴わない変更ではtasksやADRを更新しない．
- semantic decisionを伴う変更では，完了時に一度だけcurrent，tasks，ADR要否を確認する．
- 文書削除時は旧pathをまとめて検索する．
- Cloud reviewはmerge-ready candidateに対して原則1回依頼する．
- actionable feedbackはまとめて修正し，細かな再レビューを繰り返さない．

## Anti-Patterns

以下が続くとresource削減効果が下がる．

- currentへの完了履歴追記
- tasksへのPR全文コピー
- 同じparameter表の複数文書への複製
- 新規task-specific docの無制限追加
- obsolete historyを削除しない
- full CSVやlogを毎回直接読む
- Issue scopeが不明なまま重量処理を開始する
- docs-only変更でfull pytestを実行する
- decision indexやrun indexなど新しい保守台帳を安易に追加する
