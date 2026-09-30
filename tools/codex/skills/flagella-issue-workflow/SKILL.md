---
name: flagella-issue-workflow
description: Issue 単位の Codex 作業で、必要な project 規約・文書・完了手順へ段階的に案内する。
---

# Flagella Issue Workflow

## 役割

この skill は Issue 駆動作業の最小ルーターである。不変の安全・実行・完了境界は `AGENTS.md`、
詳細な Issue / PR lifecycle は `docs/codex/codex_workflow.md` を正本とする。同じ規則をここへ複製しない。

## 開始時

1. `git status --short --branch`、最新の user request、対象 Issue / PR を確認する。
2. Issue の execution target、condition 数、Mac 見積り、`execution:*` label、roadmap metadata を確認する。
3. `main` / `master` で直接作業せず、作業を planning、implementation、diagnostic、review-only、workflow、または documentation-maintenance に分類する。
4. 実装開始時は execution target、condition 数、Mac 見積り、許可された実行範囲をユーザーへ報告する。

`execution:triage`、target / label の不一致、または roadmap metadata 不備では read-only triage 以外を開始しない。
`execution:cs10` で独立conditionが2以上なら、`cs10_qualified`、condition ごとの output 分離、dry-run plan、
serial 例外のユーザー明示 IssueコメントURLを確認するまで runtime を開始しない。外部操作の承認境界は `AGENTS.md` と
`docs/codex/cs10_runbook.md` を参照する。

## 必要時の参照先

- 読む順序と Phase / run record の選択: `references/context-routing.md`
- `review_result.json`、commit、push、PR、マージ後の Issue 引継ぎ: `docs/codex/codex_workflow.md`
- 文書量・実行量の削減: `references/resource-reduction.md`
- Phase 文書の統合・移行・削除: `phase-document-maintenance` skill と `docs/codex/phase_document_policy.md`

## 完了

local review、`review_result.json`、commit、push、source-Issue-linked PR、初回完了報告の順序は
`AGENTS.md` と `docs/codex/codex_workflow.md` に従う。PR マージ後は source Issue の受入条件を確認し、
継続親Issueを子 PR だけで閉じない。残作業の新規 Issue / sub-issue は、具体案を提示してユーザー承認を得てから作成する。
