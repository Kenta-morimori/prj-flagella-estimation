# Codex Cloud connector review migration

## 目的

Codex PR reviewをrepository-managed GitHub Actions方式からChatGPT/Codex Cloud connector方式へ移行する。

## 背景

`codex-review` workflowは `OPENAI_API_KEY` を使って `openai/codex-action@v1` を実行していたが、PR #60での再実行時に `Quota exceeded. Check your plan and billing details.` で失敗した。
一方で、PR #67では `chatgpt-codex-connector[bot]` がCloud connector経由でレビューを投稿できており、自作workflowとCloud connectorの二重運用になっていた。

## 変更

- `.github/workflows/codex-review.yml` を削除する。
- Actions版専用prompt `.github/codex/prompts/review.md` を削除する。
- `AGENTS.md` のPR review方針をCloud connector前提へ更新する。
- ADR 0005をSupersededにし、ADR 0006を追加する。

## 確認

- 参照残り確認: `rg -n "codex-action|OPENAI_API_KEY|codex-review.yml|GitHub Action PR review|GitHub Actions方式|\\.github/codex/prompts/review\\.md" .github AGENTS.md docs/adr docs/codex-runs docs/planning docs/PROJECT_PLAN.md`
  - 残参照は過去ログ、Superseded ADR、新ADR、今回の移行ログのみ。
- `.github`確認: `find .github -maxdepth 3 -type f -print | sort`
  - `.github/workflows/ci.yml` のみ残存。
- `git diff --check`: pass
- `uv run ruff check .`: pass
- `uv run pytest -q`: 159 passed
