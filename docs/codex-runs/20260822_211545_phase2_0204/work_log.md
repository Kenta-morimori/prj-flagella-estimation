# Issue #204: feature-study reference run archive

- 2026-08-22 JST: `outputs/2026-08-18/174828`（36条件、約24 GB）を固定名のPhase 2 artifact pathへ移動した。
- 旧timestamp pathは相対symlinkとして復元し、既存manifestと過去の解析導線を互換維持した。
- Phase 2 / 4 currentから参照するlive reference文書を追加した。raw artifactはGit管理しない。
