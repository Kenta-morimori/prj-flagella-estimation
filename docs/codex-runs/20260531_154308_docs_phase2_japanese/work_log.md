# docs/phase2 日本語化方針の明文化

## 変更内容

- `docs/phase2/geometry_contract.md`, `docs/phase2/body_only_baseline.md`, `docs/phase2/hook_gate.md` の見出し・説明・表記を日本語中心に調整した。
- `AGENTS.md` の Language policy に、`docs/phase*/` 配下などユーザー向けプロジェクト文書は日本語で書く方針を追加した。
- 技術識別子、列名、config key、command、filename、accepted domain terms は精度維持のため原語を残す方針も明文化した。

## 検証

- `git diff --check`
- commit hook:
  - `uv run ruff format --check .`
  - `uv run ruff check .`
  - `uv run pytest -q`

## 結果

- PASS
- 対象変更 commit: `278c80c298013799de3a00e2cb13177798af6371`
