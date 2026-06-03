# docs/phase2 日本語化方針の明文化

## 変更内容

- `docs/phase2/phase2_2_geometry_contract.md`, `docs/phase2/phase2_3_body_only_baseline.md`, `docs/phase2/phase2_4_hook_gate.md` の見出し・説明・表記を日本語中心に調整した。
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
- 対象変更 commit: `54f54c453b3659a0a2207a3a4cce69a36e79000b`
