# ADR 0006: Codex Cloud Connector PR Review

## Status

Accepted

## Context

ADR 0005では、PRコメントの `@codex review` をGitHub Actionsで受け、`openai/codex-action@v1` と `OPENAI_API_KEY` でCodexレビューを実行する方式を導入した。
しかし、実運用では以下の問題が確認された。

- `OPENAI_API_KEY` のquota / billing設定に依存し、PR review workflowが `Quota exceeded` で失敗した。
- ChatGPT/Codex Cloud connectorを有効化すると、`chatgpt-codex-connector[bot]` によるPR reviewが別経路で実行され、二重レビューになった。
- repository-managed workflow、structured output schema、GitHub Review API投稿処理を保守する必要があり、運用コストが高い。

## Decision

CodexによるPR reviewは、repository-managed GitHub Actionsではなく、ChatGPT/Codex Cloud connectorへ一本化する。

- `.github/workflows/codex-review.yml` は削除する。
- GitHub Actions Secretsの `OPENAI_API_KEY` をPR review目的では使わない。
- PR reviewの起動はCloud connector側の機能に任せる。
- PR上で必要に応じて `@codex review` をコメントする。
- フィードバック対応を依頼する場合は、PRコメントで `@codex address that feedback` を使う。
- Codex reviewは補助レビューであり、task完了の正本は引き続き `docs/codex-runs/<run-id>/review_result.json` とする。

## Consequences

GitHub Actions上のAPI quota failureと二重レビューを避けられる。
また、repository側でCodex review用のstructured output schemaや投稿scriptを保守する必要がなくなる。

一方で、Cloud connectorの細かなprompt、出力schema、inline suggestion投稿ロジックはrepository内で完全には制御しない。
レビュー運用の安定性はChatGPT/Codex Cloud側の設定と権限に依存する。

## Amendment

ADR 0008 により，repository 側には Cloud connector review の実施確認だけを行う `codex-review-gate` を追加する。この gate は Codex review を実行せず，`OPENAI_API_KEY` や `openai/codex-action` も使わない。
