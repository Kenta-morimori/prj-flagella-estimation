# Issue #268: 指示の静的棚卸し

2026-10-19 時点の `main` を確認した。行数は context 量の代理指標であり、token 使用量ではない。`phase-document-maintenance/SKILL.md` は現時点で231行（Issue 作成時の179行から変化）、詳細 checklist は261行。

| 分類 | 元の問題と正本 | 対象 task・根拠 | 判断 |
| --- | --- | --- | --- |
| 恒久安全境界 | `AGENTS.md` の cs10 操作別承認、remote script 承認、秘密情報、merge 境界。`docs/codex/cs10_runbook.md` は詳細手順。 | #263 の runtime 記録と `tests/test_cs10_parallel_policy.py` が具体的な境界を要求する。 | 保持。文字列重複を削る候補から除外する。 |
| 研究・再現性 contract | pending model の `development_evaluation`、short-screen QC、replay は `AGENTS.md` と `model-development-evaluation` skill から contract へ案内する。 | Phase 2 pending model ケース。42行の skill は適用範囲が狭い。 | 保持。今回移動・削除しない。 |
| 必要時 routing | `AGENTS.md`、Issue skill、context-routing reference に phase document policy と phase skill への導線が重なる。Phase 文書の保持基準は `docs/codex/phase_document_policy.md`。 | Phase 文書整理ケース。現行 skill root 231行と詳細 checklist 261行に手順の重なりがある。 | root skill の詳細だけを短い導線にする候補を同一モデルの評価対象とする。安全規則と適用条件は保持する。 |
| 運用手順 | 初回完了報告、source Issue、targeted test が `AGENTS.md`、Issue skill、workflow、policy test に現れる。詳細の正本は `docs/codex/codex_workflow.md`。 | Issue・PR 完了ケース、docs-only ケース。`tests/test_codex_completion_reporting_policy.py` は複数文書の同一語句を要求する。 | 今回は一群ずつ評価するため変更しない。後続の候補とし、正本参照と意味の保持を検証する。 |
| 過去モデル向け補助 | 月次手順の `gpt-6-sol` 記述は #254 の初期設定判断として残る。現行設定の正本は `.codex/config.toml`。 | #268 モデル移行。 | 履歴を保持し、現行説明と例外判断のみ更新する。 |

Repository skill の自動検出先 `.agents/skills` に現存するのは `model-development-evaluation` のみ。Issue skill と Phase 文書 skill は `tools/codex/skills` にあり、`AGENTS.md` から明示参照される。Codex の利用可能 skill 一覧も前者のみを repository skill として表示した。自動検出への移動・symlink は適用過多と参照切れを招き得るため、適用・非適用ケースの評価結果を見て判断する。

PR #258 で一度欠落した `curl ... | sh` 承認と見出し参照、cs10 操作別承認、local PASS と PR / merge gate、Phase 2 QC・再現性は回帰条件とする。どれか一つでも退行する指示変更は採用しない。
