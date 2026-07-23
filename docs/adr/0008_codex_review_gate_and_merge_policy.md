# ADR 0008: Codex review gate and scoped merge policy

## Status

Accepted

## Context

`main` には branch protection / ruleset がなく，CI や `@codex review` は merge 条件として強制されていなかった。

一方，ADR 0006 により Codex review 自体は ChatGPT/Codex Cloud connector へ一本化しており，repository-managed `openai/codex-action` workflow は再導入しない方針である。

## Decision

Codex review の実行は引き続き Cloud connector に任せる。repository 側には，head SHA を含む `@codex review <head-short-sha>` が要求され，allowlisted Cloud connector が反応したことと，feedback 後は current connector-authored review threads がすべて resolved または outdated であることを検査する lightweight gate を置く。

- `.github/workflows/codex-review-gate.yml` を追加する。
- workflow は checkout せず，PR metadata / comments / reviews / reactions だけを読む。
- PR comments と PR reviews はページネーションして全件走査する。
- `OPENAI_API_KEY` や `openai/codex-action` は使わない。
- workflow は PR head SHA に commit status `codex-review-gate` を付与する。
- review request comment 自体にも現在の head SHA を含むことを要求する。
- connector response は review request comment の更新時刻以後の response だけを有効とする。
- Cloud connector の login は `chatgpt-codex-connector` / `chatgpt-codex-connector[bot]` の完全一致allowlistだけを許可する。
- Cloud connector response は，PR review の対象commit，または comment / review body / thumbs-up対象requestに含まれる head SHA で現在の head を確認する。
- Cloud connector の問題なし応答は `Final verdict PASS`，`APPROVED` review，または no-major-issues comment を成功シグナルとして扱う。
- thumbs-up reaction は review request comment の更新時刻以後の reaction だけを有効とする。
- Cloud connector が PR review だけで応答した場合に備え，`pull_request_review` event でも再評価する。
- Cloud connector が問題なしを `👍` reaction だけで返す場合があるため，`schedule` と `workflow_dispatch` でも再評価できるようにする。
- `schedule` 再評価では open PR をページネーションして走査する。
- `main` ruleset では `test` と `codex-review-gate` を required status checks にする。
- ruleset `main-required-ci-and-codex-review` は，workflow が default branch に入るまで `disabled` で作成し，merge後に `active` へ切り替える。

Codex は，小タスクに限り以下を満たす場合だけ merge してよい。

- `docs/codex-runs/<run-id>/review_result.json` が `PASS`
- CI が pass
- `codex-review-gate` が pass
- ユーザー定性評価や重大な設計判断が不要

Cloud review は merge 直前の final candidate に対して原則1回だけ要求する。指摘が出た場合は，未解決の actionable review threads を一括取得し，まとめて修正し，merge-final self-check を再実行してから全 current actionable threads を resolve する。修正不要と判断して resolve する場合は理由コメントを残す。Cloud review を commit ごとの lint 代わりに使わない。

通常 commit の hook は ADR 0007 の lightweight checks のままとし，merge-final self-check で docs/schema/Issue状態/必要テストを厚く確認する。merge-final self-check は品質を守るための手順であり，pre-commit hook を重くする理由にはしない。

重大判断やユーザー定性評価が必要な場合は，`review_result.json` を `FAIL` とし，必要な実行コマンド，出力先，評価観点，未判断点を提示して止める。

## Consequences

小さな docs / workflow / test補強 / 狭いbug fix は，CI と Codex review gate が通れば短い報告単位で進められる。

物理モデル，dataset採用条件，Phase境界，ML学習条件，動画の自然さなどの判断は自動mergeしない。

`codex-review-gate` workflow は default branch に入るまで required check として安定運用できないため，ruleset の有効化はこのADRを含むPRのmerge後に行う。
