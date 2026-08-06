# Phase Document Policy

この文書は，`docs/phase*/`に配置する文書の役割，更新条件，保持・削除基準を定める正本である．

対象はPhase 2に限定せず，Phase 3，Phase 4以降にも適用する．

## Purpose

Phase文書は，以下を両立させる．

* 現在の作業状態を短時間で把握できる
* 過去の重要判断を後続開発で再利用できる
* 現在の実装に必要なschema・contractを確認できる
* 同じ情報を複数文書へ重複させない
* 完了済みplanや一時的reportを蓄積し続けない
* Codexが必要な文書だけを段階的に読める

完全な時系列履歴は，Issue，PR，Git履歴，ADR，`review_result.json`に保持する．

通常の開発で必要な判断は，Git履歴を調べなくても`phaseX_tasks.md`から辿れる状態にする．

## Standard Structure

各Phaseは，原則として以下の構成を使う．

```text
docs/phaseX/
├── phaseX_guide.md
├── phaseX_current.md
├── phaseX_tasks.md
└── phaseX_*.md
```

`phaseX_guide.md`は，Phase固有の恒常的な規則がある場合に作成する．
規則が少ないPhaseでは省略してよい．

## `phaseX_guide.md`

Phase内で繰り返し適用する，安定した作業規則を記載する．

対象例:

* 標準CLI形式
* 入出力の基本境界
* Phase固有の長時間実行方針
* 目視レビュー観点
* モデル，dataset，trainingの解釈上の注意
* 正本となるschema，config，scriptへの導線

以下は記載しない．

* Active Issue
* 個別runの結果
* 一時的blocker
* 完了済みタスクの履歴
* task固有のacceptance criteria
* output pathや実行日時

これらは`current`，`tasks`，Issue，PRへ分離する．

## `phaseX_current.md`

現在の作業を開始するための短い入口とする．

### 保持する内容

* Phase goalの短い説明
* current baseline
* Active work
* Next queue
* current blockers
* context routing
* 現在必要な参照先

### 運用基準

* Active workは原則1件とする
* 意図的な並列作業では複数件を許容する
* Next queueは最大3件とする
* 完了したタスクは削除する
* 新しい情報の追加時に古い記述も削除する
* 目安として80行・4,000文字程度に保つ

上限よりも意味の明確さを優先するが，詳細は他文書へ分離する．

### 記載しない内容

* 完了済みIssueの詳細
* 全open Issue一覧
* 長い実験履歴
* 全数値結果
* 全実行コマンド
* output pathや実行日時
* branch履歴
* configの完全な複製
* tasksやADRの長い要約

## `phaseX_tasks.md`

後続開発で必要となる判断を，判断単位で圧縮して保持する．

Issue単位のwork logではない．
複数Issue，PR，実験が同じ判断へ収束した場合は，1つのsectionへ統合してよい．

各entryは，以下を説明する．

1. どのような問題・背景があったか
2. モデル，条件，schema，pipelineをどう変更したか
3. どの比較条件を評価したか
4. どの結果が意思決定に影響したか
5. 結果をどう解釈したか
6. 何を採用，不採用，置換，保留としたか
7. 詳細な根拠をどこから確認できるか

### 標準形式

```markdown
### <Decision ID>: <Title>

- **Status:** adopted / rejected / replaced / pending / diagnostic
- **Background:** 変更または検証が必要になった背景．
- **Change:** モデル，条件，schema，pipelineの変更内容．
- **Comparison:** 判断に必要だった比較条件．
- **Result:** 意思決定に影響した主要結果．
- **Interpretation:** 結果の考察と限界．
- **Decision:** 採用，不採用，置換，保留した内容．
- **Evidence:** Issue，PR，ADR，test，config，schema，report．
```

`Comparison`が不要な場合は省略してよい．

### Status

* `adopted`: 現在採用している
* `rejected`: 評価したが採用しなかった
* `replaced`: 過去の採用判断を別判断へ置換した
* `pending`: 評価またはユーザー判断が未完了
* `diagnostic`: 診断には有用だがcanonicalではない

### Numerical Results

残す数値:

* 採用default
* pass / failを分けた主要境界
* 採否を決めた代表比較値
* 論文値との差を説明する代表値
* 後続実装で誤解されやすい重要値

原則として残さない数値:

* 全seedの結果
* 全sweep cell
* 全時刻
* 全output path
* 判断に影響しなかった診断値
* 同じ結論を示す重複値

詳細値はIssue，PR，ADR，config，test，report，Git履歴から辿れるようにする．

### Excluded Content

tasksへ通常記載しないもの:

* branch名
* 完全なacceptance criteria
* 全実行コマンド
* 実装ファイル一覧
* output directory一覧
* 実行日時
* PR descriptionの複製
* 時系列の作業日誌

## Task-Specific Documents

個別文書は，現在の実装・運用で直接参照する場合に保持する．

### 保持対象

* schema definition
* output-column definition
* feature definition
* validation contract
* parameter mapping
* dataset version registry
* Phase handoff contract
* active experiment / validation procedure
* 後続解析で再利用する大規模report
* tasksへ圧縮すると実装情報が欠損する表や対応関係

### 統合・削除候補

* 完了済みplan
* temporary design memo
* 過去候補の列挙
* PR本文と重複する実装記録
* 実行コマンドとoutput pathが中心の文書
* tasksとADRだけで判断を再現できる文書
* 現在のcode，schema，pipelineから参照されない文書

ファイル名だけで保持・削除を判断しない．
本文，参照元，config，test，schema，ADRを確認する．

## ADR Boundary

ADRは，詳細な設計理由を独立して残す必要がある判断に使用する．

主な対象:

* 物理解釈の変更
* reference modelへのproject-specific extension
* architecture変更
* schemaまたは共通contract変更
* dataset採択policy変更
* Phase boundary変更
* ML training / evaluation policy変更
* 既存の重要設計の置換

`phaseX_tasks.md`には結論と代表根拠だけを記載し，ADR本文を複製しない．

## When To Update Tasks

以下を採用，変更，撤回，置換した場合はtasksを更新する．

* physical model
* geometry interpretation
* standard execution condition
* numerical stabilization policy
* validation / quality gate
* output schema
* dataset versionまたは採択範囲
* Phase handoff
* pipeline boundary
* data mixing policy
* ML training / evaluation policy
* 重要なfailure modeの解釈

通常，以下では更新しない．

* typoやformatting
* semantic effectのない内部refactor
* contractを変えないtest cleanup
* 解釈を変えないlogging変更
* supported interfaceを変えないcommand例の修正

## Completion Lifecycle

semantic decisionを伴うタスク完了時は，以下を行う．

1. `phaseX_current.md`から完了タスクを外す
2. baseline，Next queue，blockerを更新する
3. `phaseX_tasks.md`へ判断を追加・更新する
4. ADRの作成・更新要否を判断する
5. task-specific docを分類する
6. obsolete historyを統合・削除する
7. 移動・削除した文書への参照を更新する
8. `review_result.json`へ文書更新内容を記録する

task-specific docは以下に分類する．

* `live contract`
* `active validation`
* `reusable report`
* `obsolete history`

## Deletion Gate

文書を削除する前に，以下をすべて確認する．

* 背景，変更，主要結果，考察，採否がtasksまたはADRへ移行済み
* live contract情報がconfig，schema，test，または別文書に保持されている
* 詳細な証跡をIssue，PR，ADR，report，Git履歴から確認できる
* 削除対象への参照が更新されている
* currentとtasksの採用状態が一致している

Git履歴に残ることだけを理由に，通常開発で必要な判断を省略しない．

## Source Of Truth

| 情報                        | 正本                       |
| ------------------------- | ------------------------ |
| 現在の作業状態                   | `phaseX_current.md`      |
| 採択判断と考察                   | `phaseX_tasks.md`        |
| Phase固有の恒常規則              | `phaseX_guide.md`        |
| 詳細な設計理由                   | ADR                      |
| parameter値                | config                   |
| machine-readable field    | schema                   |
| gate実装                    | codeとtest                |
| scopeとacceptance criteria | Issue / PR               |
| 完全な作業履歴                   | Issue / PR / Git history |
| Codex実行結果                 | `review_result.json`     |
| raw diagnostic result     | generated output         |

同じ情報を複数箇所へ記載する場合は，正本と要約の関係を明確にする．

## Cross-Phase Application

このpolicyは全Phaseに適用する．

Phase固有の差分は`phaseX_guide.md`またはlive contractへ記載し，共通規則を複製しない．

Phaseを跨ぐ判断は，以下のいずれかへ記録する．

* 上流Phaseのhandoff contract
* 下流Phaseのinput contract
* ADR
* 両Phaseのtasksへの短い相互参照

## Verification

最低限，以下を確認する．

```bash
git diff --check
```

文書を削除・名称変更した場合は，旧pathと旧filenameを検索する．

```bash
rg -n "<deleted-path>|<deleted-file-name>" \
  AGENTS.md README.md docs tools scripts conf schemas tests \
  -g'!docs/codex-runs/**'
```

追加確認:

* 参照pathが存在する
* current，tasks，ADR，config，schema，testの用語が一致する
* deprecated terminologyが正本として残っていない
* docs-onlyでない場合は対応するtargeted testを実行する

docs-only，planning-only，workflow-only変更ではfull pytestを必須としない．

## Anti-Patterns

以下を避ける．

* currentへの完了履歴追記
* tasksへのIssue / PR全文コピー
* tasksへの全sweep結果貼り付け
* 同じparameter表の複数文書への複製
* 一時的planを恒久contractとして残す
* active validation完了後に保持要否を再評価しない
* 削除文書への参照を残す
* decision indexやhistory文書を追加して保守対象を増やす
* Git履歴だけを通常の判断記録として扱う
