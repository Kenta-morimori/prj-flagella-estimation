---
name: flagella-issue-workflow
description: Use for issue-driven Codex work in this repository when handling GitHub issues, task planning, implementation, simulation/test execution, result analysis, review_result logging, commits, pushes, PR creation, or when reducing context/resource usage by routing only the necessary project documents.
---

# Flagella Issue Workflow

## Overview

この skill は、`prj-flagella-estimation` の issue 単位作業を、必要な文書だけを読んで進めるための workflow である。
目的は、`AGENTS.md` の完了条件を守りながら、方針検討・実装・実行・分析・報告の各段階で context と実行コストを使いすぎないことである。

## Start

1. 最初に `git status --short --branch` と現在 branch を確認する。
2. issue / PR / user request の対象を特定する。
3. `main` / `master` 上で直接作業しない。新規 issue 対応は最新 `main` から feature branch を作る。
4. まず作業タイプを分類する。

作業タイプ:

- `planning`: 方針検討、タスク分解、issue draft、実装前整理。
- `implementation`: コード、設定、テスト、ドキュメントの変更。
- `diagnostic`: Phase 2 simulation sweep、失敗条件保存、原因切り分け。
- `review-only`: 差分レビュー、CI確認、既存結果の整理。
- `workflow`: Codex 運用、skill、review_result、PR作成方針の変更。

## Context Routing

常に全部読まず、作業タイプで読む文書を絞る。

- どの文書を読むか迷う場合は `references/context-routing.md` を読む。
- review_result、commit、push、PR 完了条件が必要な場合は `references/completion-policy.md` を読む。
- リソース削減方針、軽量/重量タスクの分け方、今後の issue 化候補が必要な場合は `references/resource-reduction.md` を読む。

基本ルール:

- `AGENTS.md` は repo policy の正本として扱う。
- `docs/PROJECT_PLAN.md` は phase 全体の現在地を確認するときだけ読む。
- `docs/phase2/phase2_tasks.md` は採択済み Phase 2 task の status / acceptance criteria を確認するときだけ読む。
- 個別 task の詳細は、該当する `docs/phase*/phase*_*.md` だけ読む。
- 過去 run log は、該当 issue / branch / phase に直接関係するものだけ読む。

## Workflow

1. **方針検討**
   - issue 本文、関連 task、直近の該当 docs だけで scope を決める。
   - 実装に入る前に、完了条件、実行すべき test/simulation、目視レビュー要否を分ける。
   - 方針検討だけの依頼では、重い simulation や full pytest は走らせない。

2. **実装**
   - 変更は scope 内に限定する。
   - `scripts/` は orchestration、再利用ロジックは `src/`、Codex workflow 補助は `tools/codex/` に置く。
   - Phase 2 の物理モデル変更では、参照論文モデル、現行実装、数値安定化、project-specific extension を区別して記録する。

3. **実行**
   - 最小の targeted test から実行する。
   - Phase 2 の simulation は、まず短時間または代表条件で確認し、必要があるときだけ sweep / 長時間実行へ進む。
   - `time.dt_star` は default config を変えず、必要な run で CLI override する。

4. **結果分析**
   - PASS / FAIL だけでなく、残った blocking issue と次に切り分けるべき原因を記録する。
   - Phase 2 の collapse、fly-away、hook drift、no_bundle などは、失敗でも再現条件が有用なら diagnostic progress として扱う。

5. **報告・完了**
   - `docs/codex-runs/<run-id>/review_result.json` を作成する。
   - PASS 完了なら commit / push / PR まで行う。
   - FAIL 診断でも有用な場合は、diagnostic / wip / docs commit として保存できる。ただし完了とは言わない。

## Resource Discipline

- `SKILL.md` に長い project history を入れない。長い情報は reference へ分離する。
- 1回の issue で「方針検討」と「実装・実行」を混ぜすぎない。必要なら PR を分ける。
- full test / long simulation / video render は、acceptance criteria に必要な場合だけ実行する。
- 既存 docs の探索は `rg` で候補を絞ってから読む。
- 最終報告では、読んだもの、変えたもの、走らせたもの、走らせなかったものを短く明示する。
