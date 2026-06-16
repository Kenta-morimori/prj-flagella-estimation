# Codex Run Log

## Task

Issue #64「CodexによるPRレビュー」に対応し、PRコメントの `@codex review` でCodex GitHub Actionを起動するworkflowを追加する。
Copilot風の運用に寄せるため、可能な指摘はGitHub Review APIのinline suggestionとして投稿し、inline化できない指摘はReview本文にfallbackする。

## Scope

- `.github/workflows/codex-review.yml` を追加する。
- `.github/codex/prompts/review.md` を追加する。
- `AGENTS.md` にCodex GitHub Action PRレビュー方針を追記する。
- Codex workflow変更としてADRを作成する。

## Implementation Notes

- `issue_comment.created` で起動し、対象がPRかつコメント本文に `@codex review` がある場合だけCodex jobを実行する。
- `openai/codex-action@v1` は `sandbox: read-only` と `safety-strategy: drop-sudo` で実行する。
- Codexには構造化JSONを返させ、GitHubへの投稿は固定の `actions/github-script` で行う。
- `path` と `line` がPR diff上の変更行に対応できる場合だけinline commentとして投稿する。
- `suggested_code` はGitHub Markdownの `suggestion` blockへ変換する。
- 行対応できない指摘、JSONの軽微な崩れ、上限超過はReview本文のfallback notesへ回す。
- レビュー本文は日本語、Review本文の見出しと `Final verdict` の `PASS` / `FAIL` は固定形式にする。

## Checks

- `git diff --check`
- `ruby -e 'require "yaml"; YAML.load_file(".github/workflows/codex-review.yml"); puts "yaml ok"'`
- `command -v actionlint` は未検出のため、ローカルではactionlint未実行。
- `uv run ruff format --check .`
- `uv run ruff check .`
- `uv run pytest -q`
- `rg -n "OPENAI_API_KEY\\s*[:=]|sk-[A-Za-z0-9]|github_pat_|ghp_" .github AGENTS.md docs/adr/0005_codex_pr_review_action.md`

## Review

ユーザー目視レビューは不要。
この変更はPhase 2物理モデルや動画出力を変更しないため、Phase 2 visual review policyの対象外である。
