# ADR 0005: Codex PR Review Action

## Status

Accepted

## Context

Issue #64では、PRコメントで `@codex review` が投稿された時だけCodex CLIによるレビュー補助を実行したい。
従来のCopilot CLI補助に近い運用として、全体コメントだけでなく、可能な箇所ではGitHubのinline suggestionとして具体的な修正案を提示したい。

一方で、GitHub Actions上でOpenAI API keyを扱うため、PR本文・差分・コメントからのprompt injection、外部ユーザーによるAPI使用量悪用、secret漏洩を避ける必要がある。
また、`openai/codex-action@v1` はCodexの最終出力を返すActionであり、PR review commentやcommit suggestionを直接投稿する機能は持たない。

## Decision

Codex PRレビューは以下の構成で導入する。

- 起動条件は `issue_comment.created` とし、PRコメントに `@codex review` が含まれる場合だけ実行する。
- 実行許可は `openai/codex-action` のwrite access checkに任せ、外部ユーザー全員を許可しない。
- Codexは `sandbox: read-only` と `safety-strategy: drop-sudo` で実行する。
- Codexにはレビュー結果を構造化JSONとして出力させる。
- GitHubへの投稿はCodexに任せず、workflow内の固定 `actions/github-script` で検証してから `pulls.createReview` を呼ぶ。
- `path` と `line` がPR diff上の変更行に対応できる指摘だけinline review commentとして投稿する。
- `suggested_code` がある場合だけGitHub Markdownの `suggestion` blockへ変換する。
- inline化できない指摘はReview本文のfallback notesへ回す。
- レビュー本文は日本語、見出しと `PASS` / `FAIL` は固定形式にする。

## Consequences

この方式では、Copilotに近いinline suggestion体験を提供しつつ、CodexにGitHub write操作を直接任せない。
行番号やdiff対応が不確かな指摘は本文にfallbackするため、レビュー投稿の失敗や誤った行へのsuggestionを減らせる。

一方で、GitHub Review APIに渡せるのはdiff上の有効な行だけである。
そのため、ファイル全体の設計指摘、テスト不足、複数行・複数ファイルの修正提案はinline suggestionではなくReview本文に集約する。

このPRレビューActionの `PASS` / `FAIL` はPR補助レビューの判定であり、リポジトリ作業完了の正本である `docs/codex-runs/<run-id>/review_result.json` とは別物として扱う。
