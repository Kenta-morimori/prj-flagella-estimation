# Phase 2 task reorganization after bundling review

## Summary

Phase 2.7 の複数べん毛束化検証を完了扱いにし、次の中心タスクを Issue #65 と Issue #69 に整理した。

## Inputs

* User instruction: トルク伝搬の拡張モデルにおける束化検証は完了。続くタスクは #65, #69 が中心。
* Related issues:
  * `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/10`
  * `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/53`
  * `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/65`
  * `https://github.com/Kenta-morimori/prj-flagella-estimation/issues/69`

## Changes

* `docs/phase2/phase2_current.md`
  * Phase 2.6 / 2.7 の完了状態を反映。
  * 現在の主対象を Issue #65 と Issue #69 に更新。
  * hook length fail を #65 の前提リスクとして明記。
* `docs/PROJECT_PLAN.md`
  * Phase 2.7 を完了に更新。
  * Phase 2.8 を RUN固定状態の本数差評価へ更新。
  * Phase 2.9 を Tumble状態の段階実装へ更新。
* `docs/phase2/phase2_tasks.md`
  * P2-6-009 と P2-7-006 を complete に更新。
  * P2-8-008 を Issue #65 の RUN本数差評価として整理。
  * P2-9-009 を完了済みsupport taskへ移動。
  * P2-9-010 として Issue #69 の段階Tumble実装を追加。
* GitHub issues
  * Issue #65 の本文を RUN固定本数差評価として更新。
  * Issue #69 の本文を段階Tumble実装として更新。
  * Issue #10 と #53 に、#65/#69 へ続く整理コメントを追加。

## Checks

* `git diff --check`
* `python -m json.tool docs/codex-runs/20260617_162732_phase2_task_reorg_after_bundling/review_result.json`
* `rg -n "#65|#69|P2-7|P2-8|P2-9|Tumble|RUN|hook_wrapped_axis_aligned" docs/phase2 docs/PROJECT_PLAN.md`
* `wc -l docs/phase2/phase2_current.md`
* GitHub API read checks for Issue #65 and Issue #69 bodies

## Notes

This is a documentation and task organization change. Physics model, simulation code, config behavior, and output formats were not changed.
