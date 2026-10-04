# Resource Reduction

この reference は、context と実行量を抑えるための入口である。文書の配置・保持基準は
`docs/codex/phase_document_policy.md`、repository 共通規約は `AGENTS.md` を正本とする。

- `rg -n` で長い Markdown、log、CSV、generated output の対象を絞る。
- `phaseX_current.md`、関連する `phaseX_tasks.md`、live contract / config / schema / test の順に読み、
  前段で足りる場合は ADR、Git history、raw output を読まない。
- 過去 run は `review_result.json` を先に読み、必要な場合だけ `work_log.md` を読む。
- docs-only / workflow-only 変更で full pytest を既定にせず、最小の relevant check を実行する。
- 長時間 simulation、sweep、training、render は、受入条件とユーザー明示依頼がある場合だけ実行する。
- Phase 文書の統合・削除は `phase-document-maintenance` skill を使い、旧 path の stale reference を検索する。

判断根拠、採択・不採択・保留、現行 schema / contract / config / test、再現性に必要な artifact は削減対象にしない。
