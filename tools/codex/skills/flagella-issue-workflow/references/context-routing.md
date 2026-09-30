# Context Routing

常に `AGENTS.md`、最新 user request、対象 Issue / PR、`git status --short --branch` を確認する。
以後は次の目的に応じて必要な正本だけを読む。

| 確認したいこと | 読む正本 |
| --- | --- |
| repository の安全、実行、完了境界 | `AGENTS.md` |
| Phase の現在地と次に読む文書 | `docs/phaseX/phaseX_current.md` |
| Phase 固有の運用規則 | `docs/phaseX/phaseX_guide.md` |
| 過去の採択・不採択・保留 | `docs/phaseX/phaseX_tasks.md` を `rg -n` で検索 |
| 現行の parameter、field、gate | config、schema、code、test、active validation 文書 |
| 重要な設計理由 | 関連 ADR |
| Issue / PR の scope、受入条件、review | source Issue、target PR、review comment |
| 過去 run の状態 | `review_result.json`、必要時のみ `work_log.md` |
| 文書の統合・削除 | `phase-document-maintenance` skill と `docs/codex/phase_document_policy.md` |
| commit、push、PR、merge、引継ぎ | `docs/codex/codex_workflow.md` |

Issue / PR、run log、generated output を Phase 文書へ全文複製しない。docs と machine-readable source が矛盾する場合は、
無断で一方へ合わせず報告する。
