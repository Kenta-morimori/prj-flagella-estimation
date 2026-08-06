# Phase Document Consolidation Checklist

このreferenceは，`phase-document-maintenance` skillでPhase文書を統合・削除するときに使用する詳細チェックリストである．

文書構成と保持基準の正本は`docs/codex/phase_document_policy.md`とする．

## 1. Scope

- [ ] 対象Phaseを特定した
- [ ] 対象Issue / PRを特定した
- [ ] auditのみか，実際に変更するかを確認した
- [ ] 文書削除が許可されているか確認した
- [ ] semanticな仕様変更を含むか確認した
- [ ] scope外のPhase，code，config，datasetを明記した

## 2. Inventory

対象Phaseの文書を以下の形式で整理する．

```markdown
| File | Classification | Decision-bearing content | Migration target | Evidence | Action |
|---|---|---|---|---|---|
| `phaseX_...md` | live contract | column定義 | keep | schema / test | retain |
| `phaseX_...md` | obsolete history | 採択理由と代表結果 | tasks / ADR | Issue / PR | migrate then delete |
```

Classification:

- `current`
- `tasks`
- `guide`
- `live contract`
- `active validation`
- `reusable report`
- `obsolete history`
- `unclear`

Action:

- `retain`
- `update`
- `migrate then delete`
- `defer`

## 3. Per-File Review

各ファイルについて確認する．

### Purpose

- [ ] この文書の元の目的を確認した
- [ ] 現在もその目的が有効か確認した
- [ ] current / tasks / guide / contractのどれに属するか判断した

### Decision-Bearing Information

- [ ] Background: なぜ変更・検証が必要だったか
- [ ] Change: model，condition，schema，pipelineをどう変えたか
- [ ] Comparison: 何と何を比較したか
- [ ] Result: どの結果が判断へ影響したか
- [ ] Interpretation: どう考察され，どの限界があるか
- [ ] Decision: 何を採用，不採用，置換，保留したか
- [ ] Evidence: Issue，PR，ADR，test，config，schema，reportを辿れるか

### Live Contract Check

- [ ] schemaまたはoutput fieldを定義していないか
- [ ] featureやmetricの意味を定義していないか
- [ ] validation gateや判定手順を定義していないか
- [ ] parameterと論文記述の対応を定義していないか
- [ ] dataset versionまたはPhase handoffを定義していないか
- [ ] 現在進行中のvalidation procedureではないか
- [ ] 後続解析で再利用する詳細表を含んでいないか

いずれかがYesの場合，tasksへの要約だけで削除してよいか再検討する．

### Duplication Check

- [ ] 同じ内容がcurrentへ重複していないか
- [ ] 同じ内容がtasksへ重複していないか
- [ ] ADR本文を複製していないか
- [ ] config値を文書内へ重複記載していないか
- [ ] schema field一覧を複製していないか
- [ ] PR descriptionまたはacceptance criteriaを複製していないか

## 4. Tasks Migration

判断情報を`phaseX_tasks.md`へ移す場合に確認する．

- [ ] Decision IDとtitleが判断単位になっている
- [ ] Issue単位の作業日誌になっていない
- [ ] 複数Issue / PRが同じ判断なら統合した
- [ ] `Status`が正しい
- [ ] `Background`がある
- [ ] `Change`がある
- [ ] 必要な場合は`Comparison`がある
- [ ] `Result`が意思決定に関係する内容へ限定されている
- [ ] `Interpretation`に考察と限界がある
- [ ] `Decision`が採用・不採用・置換・保留を明示する
- [ ] `Evidence`から詳細へ到達できる

### Numeric Compression

残す候補:

- [ ] 採用default
- [ ] pass / failを分けた主要境界
- [ ] 採否を決めた代表比較値
- [ ] 論文値との差を説明する代表値
- [ ] 後続実装で誤解されやすい値

削除候補:

- [ ] 全seed結果
- [ ] 全sweep cell
- [ ] 全時刻
- [ ] output path一覧
- [ ] 判断に影響しなかった診断値
- [ ] 同じ結論を示す重複値

数値を削除しても，詳細なevidenceを辿れることを確認する．

## 5. Current Update

- [ ] Phase goalが短く記載されている
- [ ] current baselineが明確
- [ ] Active workは原則1件
- [ ] 並列作業の場合のみ複数件
- [ ] Next queueは最大3件
- [ ] blockersが現在の依存関係だけを示す
- [ ] context routingが短い
- [ ] 完了済みIssueの詳細を削除した
- [ ] 全数値表を削除した
- [ ] 全実行コマンドを削除した
- [ ] output pathと実行日時を削除した
- [ ] configの完全な複製を削除した
- [ ] 目安として80行・4,000文字程度に収まる

## 6. Retention Decision

### Retain

以下のいずれかに該当する場合は原則維持する．

- [ ] live schema / contract
- [ ] output-column / feature definition
- [ ] parameter mapping
- [ ] dataset registry
- [ ] Phase handoff
- [ ] active validation
- [ ] reusable large report
- [ ] tasksへ圧縮すると実装情報が欠損する

### Delete Candidate

以下をすべて満たす場合に削除候補とする．

- [ ] 完了済みplan，temporary note，過去比較，実行報告である
- [ ] decision-bearing informationをtasksまたはADRへ移行した
- [ ] live contract情報が別の正本に残る
- [ ] detailed evidenceをIssue / PR / ADR / report / Git historyから辿れる
- [ ] 他文書からの参照を更新できる
- [ ] current / tasks / ADRの採用状態が一致する

### Defer

以下の場合は削除を保留する．

- [ ] active validationが未完了
- [ ] user decisionがpending
- [ ] schema / contractかどうか不明
- [ ] configやtestだけでは意味が再現できない
- [ ] Issue / PRの結論が不明確
- [ ] tasksへ移す情報の根拠が不足する

保留理由を明記する．

## 7. ADR Boundary

- [ ] 物理解釈の変更ではないか
- [ ] reference modelへのproject-specific extensionではないか
- [ ] architecture変更ではないか
- [ ] schema / shared contract変更ではないか
- [ ] dataset policy変更ではないか
- [ ] Phase boundary変更ではないか
- [ ] ML training / evaluation policy変更ではないか
- [ ] 既存の重要設計を置換していないか

該当する場合，tasksだけで済ませずADRの作成・更新要否を判断する．

## 8. Reference Update

削除・移動・名称変更したファイルについて検索する．

```bash
rg -n "<deleted-path>|<deleted-file-name>" \
  AGENTS.md README.md docs tools scripts conf schemas tests
```

必要に応じて判断語も検索する．

```bash
rg -n "<Issue-number>|<Decision-ID>|<model-name>|<config-key>|<dataset-version>" \
  docs/phaseX docs/adr conf schemas tests
```

確認項目:

- [ ] Markdown linkを更新した
- [ ] plain textの旧filenameを更新した
- [ ] context-routingの参照を更新した
- [ ] scripts / config READMEの参照を更新した
- [ ] deprecated terminologyが正本として残っていない
- [ ] 削除後の参照先がtasks / ADR / live contract / Issue / PRのいずれかである

## 9. Cross-Source Consistency

- [ ] currentとtasksの採用状態が一致する
- [ ] tasksとADRの結論が一致する
- [ ] tasksのparameter値がconfigと一致する
- [ ] tasksのfield名がschemaと一致する
- [ ] tasksのgate値がcode / testと一致する
- [ ] dataset versionとregistryが一致する
- [ ] Phase handoffのinput / outputが上下流で一致する
- [ ] pending decisionをadoptedとして記載していない

## 10. Verification

最低限:

```bash
git diff --check
```

追加確認:

- [ ] 新規skillのYAML frontmatterが有効
- [ ] skillの`name`とdirectory名が一致
- [ ] 記載されたpathが存在する
- [ ] policyとskillが矛盾しない
- [ ] 削除対象への参照検索が0件，または意図した履歴参照だけである
- [ ] docs-only変更でない場合はtargeted testを実行した
- [ ] 実行しなかったcheckと理由を記録した

## 11. Review Result

`review_result.json`へ以下を記録する．

- [ ] 更新したcurrent
- [ ] 更新したtasks
- [ ] 統合したdecision
- [ ] 維持したlive documents
- [ ] 削除したobsolete documents
- [ ] 削除を保留したdocumentsと理由
- [ ] stale-reference check
- [ ] semantic changeの有無
- [ ] tests / checks
- [ ] remaining issues

Issue作業全体の完了条件は，`flagella-issue-workflow/references/completion-policy.md`に従う．
