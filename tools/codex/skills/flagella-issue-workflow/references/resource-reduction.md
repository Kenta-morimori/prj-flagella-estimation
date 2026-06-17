# Resource Reduction

この reference は、Codex の context 使用量、tool 実行量、simulation/test 実行コストを下げるための方針である。

## Does This Skill Reduce Resources?

削減できる可能性は高い。

理由:

- 毎回読む文書を `task type` と `phase` で routing できる。
- `SKILL.md` 本体を短くし、必要時だけ reference を読む progressive disclosure にできる。
- 方針検討だけの依頼で、重い simulation や full pytest を走らせる判断を避けやすい。
- `review_result.json` と `work_log.md` の読み方を固定し、過去ログ探索を `review_result` 優先にできる。

限界:

- Skill は自動で source of truth を更新しない。`AGENTS.md`、`docs/PROJECT_PLAN.md`、`docs/phase2/*` が肥大化し続けると効果は下がる。
- 実際に効果を出すには、issue 本文に task type、対象 phase、期待する完了範囲を書く必要がある。

## Lightweight vs Heavyweight Work

軽量扱い:

- issue 整理、方針検討、task decomposition。
- docs/planning への proposal。
- review_result / work_log の整理。
- 小さな docs 変更。
- targeted unit test だけで十分な変更。

重量扱い:

- Phase 2 simulation sweep。
- `duration_s>=0.5` かつ `time.dt_star=1.0e-4` の長時間実行。
- 3D/2D動画 render。
- full pytest。
- 物理モデルや出力 format の変更。

運用:

- 軽量 task では、まず `git diff --check`、JSON/YAML validation、該当 docs の整合だけでよいか判断する。
- 重量 task では、短時間 representative、targeted tests、sweep、full tests の順で段階実行する。
- 既存PRで継続すべき Phase 2 diagnostic は、原則として同じPRブランチへ統合する。誤って別branchに進んだ場合は、差分を本来のPRブランチへ移し、重複PRは閉じる方針にする。

## Issue Template Suggestions

今後の issue には、以下を入れると resource 使用量を下げやすい。

```markdown
## Task type
- [ ] planning
- [ ] implementation
- [ ] diagnostic
- [ ] review-only
- [ ] workflow

## Scope
- Target phase:
- Target files/modules:
- Out of scope:

## Required context
- Must read:
- Read only if needed:

## Execution budget
- Tests:
- Simulation:
- Video render:
- Full pytest required: yes/no

## Completion target
- PASS completion required: yes/no
- Diagnostic FAIL can be committed: yes/no
- User visual review required: yes/no/unknown
```

## Other Reduction Ideas

1. `docs/codex-notes/phase2_current_index.md` を作る。
   - Phase 2 の current representative、最新 issue、最新 PR、読むべき doc を1ページにまとめる。
   - 長い `docs/PROJECT_PLAN.md` と `docs/phase2/phase2_tasks.md` を毎回読む頻度を下げる。

2. `docs/codex-runs/RUN_INDEX.md` を作る。
   - run-id、issue、branch、status、commit、PR、next action だけを表にする。
   - 過去 run の `review_result.json` 全探索を減らす。

3. issue を `planning PR` と `implementation PR` に分ける。
   - 方針が不確かな task では、まず docs/planning と review_result だけで PASS する軽量 PR を作る。
   - 実装 PR は acceptance criteria と execution budget が固まってから始める。

4. simulation sweep の manifest summary を標準化する。
   - sweep summary の主要列と代表条件を固定し、毎回 CSV 全体を読まなくてよいようにする。

5. 「full pytest が必要な条件」を明文化する。
   - docs-only / workflow-only では full pytest を必須にしない。
   - physics / output format / shared library 変更では full pytest または未実行理由を必須にする。

6. user visual review issue を分離する。
   - 動画生成と目視判定が必要な場合、生成 PR と review result 更新 PR を分けると無駄な再実行を減らせる。
