# ADR 0022: Codex のPRマージ後Issue引継ぎ

## Status

Accepted

## Context

source Issue を参照する PR が merge されても、PR の変更範囲、Issue の受入条件、
`review_result.json`、GitHub 上の Issue 状態は自動的に一致するとは限らない。特に継続親Issueは、
一つの child PR が merge されたことだけで完了と誤認し得る。また、残作業を即座に新規 Issue 化すると、
未承認の scope 拡張や不要な dependency を増やす。

## Decision

- repository の不変な安全・実行・完了境界は `AGENTS.md`、詳細な Issue / PR lifecycle は
  `docs/codex/codex_workflow.md`、Issue skill は必要な正本へ案内する最小ルーターとする。
- PR merge 後は、merged PR、source Issue、受入条件、relevant check、local `review_result.json` を照合する。
  すべての受入条件が満たされる場合だけ、必要な Issue checkbox・状態・完了記録を更新する。
- 継続親Issueは、child PR / child Issue の merge / close だけでは閉じず、親自身の受入条件と継続目的を確認する。
- 残作業は範囲、受入条件、`Heavy/runtime execution target` を含む具体的な後続 task として提示する。
  新規 Issue / sub-issue の作成はユーザー承認後に限る。

## Consequences

- merge 後の完了報告には、source Issue の状態と未達受入条件の有無を含める。
- GitHub 更新権限がない task では、更新内容を報告して外部状態を更新済みと主張しない。
- `AGENTS.md` と Issue skill の重複を減らすが、cs10 の操作境界と再利用契約は削除しない。
