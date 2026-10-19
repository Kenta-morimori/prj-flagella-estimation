---
name: phase-document-maintenance
description: Use when auditing, restructuring, consolidating, migrating, or deleting docs/phase*/ documents, including phase current files, tasks decision records, guides, contracts, validation documents, plans, and reports.
---

# Phase Document Maintenance

この skill は `docs/phase*/` の判断と参照を保ちながら文書を整理する。文書構成と保持基準の正本は `docs/codex/phase_document_policy.md`。詳細な inventory、判断抽出、移行、削除、検証手順は `references/consolidation-checklist.md` を必要時に読む。

## 適用と境界

`phaseX_current.md`、`phaseX_tasks.md`、guide、contract、validation、plan、report の監査・統合・移行・削除に使う。単純な typo やリンク修正には全面適用しない。

開始時に `AGENTS.md`、依頼と Issue、phase document policy、対象 Phase の current と tasks の関連箇所を確認する。長い文書は `rg -n` で先に絞る。対象 Phase、変更の有無、削除可否、semantic change の有無を明確にする。文書整理だけなら model、parameter、dataset 採択、schema、training policy を変更しない。

## 手順

1. 対象文書を inventory し、live contract・active validation・判断記録・古い履歴を区別する。
2. 統合候補から採択・却下・保留の判断と根拠を抽出し、`phaseX_tasks.md` に必要な判断だけを残す。`phaseX_current.md` は現在地と次の導線に限る。
3. live contract と後続解析に必要な文書を保持する。判断・参照先・復元可能性を確認してから古い文書を削除する。不明な文書は保留する。
4. 移動・削除後に参照元を更新し、current・tasks・ADR・config・schema・test の整合を確認する。
5. `git diff --check`、変更範囲に応じた参照切れ検索と targeted check を行う。docs-only に full pytest を一律要求しない。
6. 完了手順と `review_result.json`、commit、push、PR は `flagella-issue-workflow` と `docs/codex/codex_workflow.md` に従う。

情報移行前の削除、Git 履歴だけを根拠にした省略、pending 判断の採択扱い、根拠のない model・dataset 仕様変更はしない。
