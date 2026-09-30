# AGENTS.md

このファイルは、Codex が常に守る最小のリポジトリ規約である。詳細な Issue / PR lifecycle は
`docs/codex/codex_workflow.md`、Issue 単位の作業導線は
`tools/codex/skills/flagella-issue-workflow/` を正本とする。

## プロジェクト概要

このリポジトリは、遊泳顕微鏡動画から細菌のべん毛本数を推定するパイプラインを開発する。

1. Phase 1: repository、CLI、config、logging、再現性の基盤
2. Phase 2: 3D physical simulation model、数値・物理検証、長時間安定性、canonical model の凍結
3. Phase 3: 2D projection、pseudo-microscopy、observability、detection、細胞単位 clip 生成、dataset の凍結
4. Phase 4: flagella-count model の学習と評価
5. Phase 5+: 予測可視化と実データ解析支援

## コンテキストの読み分け

必要なものだけを、次の順で読む。

1. `AGENTS.md`
2. ユーザーの最新依頼と対象 Issue / PR
   - Issue の `Heavy/runtime execution target` と `execution:*` label を確認する。
   - 不在、不一致、または `execution:triage` なら、triage 完了まで read-only 調査に限る。
3. `docs/phaseX/phaseX_current.md`
4. 必要な場合だけ `docs/phaseX/phaseX_guide.md`
5. 過去判断が必要な場合だけ `docs/phaseX/phaseX_tasks.md` の関連箇所
6. 現行 schema、contract、config、test、active validation 文書
7. ADR
8. Issue / PR 履歴、Git 履歴、`docs/codex-runs/*/review_result.json`

長い Markdown、log、CSV、generated output を開く前に `rg -n` で対象を絞る。過去 Codex run は
`review_result.json` を `work_log.md` より先に読む。`outputs/` 配下の大きなファイルは compact summary と
manifest で不足する場合だけ読む。Phase 2 diagnostics では `run_summary.json` を先に読み、通常の分析で
`step_summary.csv` 全体を読み込まない。

ユーザーが #252 の月次棚卸し、既定 model / reasoning effort の変更提案、または AGENTS / skill の大規模再編を
明示した場合だけ `docs/codex/monthly_harness_review.md` を読む。Issue 作成用 task は調査・設計・Issue 操作・
引継ぎまでとし、repository 変更と PR 作成は実装用 task で行う。通常の実装依頼では月次手順を読まない。

## 言語

* ユーザーとは、明示的な指定がない限り日本語でやり取りする。
* ユーザー向け project 文書は既定で日本語にする。
* 技術識別子は、翻訳すると精度が下がる場合に原文を保つ。

## リポジトリ共通規約

* `main` / `master` で直接作業しない。変更前に branch と `git status` を確認する。
* 依頼範囲内に変更を限定し、必要のない大規模 refactor や dependency 追加をしない。secret、token、credential、private data、生成した認証ファイルを commit しない。
* target branch は task / Issue 指定を優先し、なければ default branch とする。実装前に Issue execution target、独立 condition 数、Mac wall time 見積り、許可された実行範囲を短く報告する。
* `execution:cs10` では、独立 condition が2以上なら `cs10_qualified` parallel-job YAML、condition ごとの output 分離、dry-run plan を確認するまで実行しない。serial 例外には Issue runbook の具体的理由と、開始前のユーザー明示 Issue コメント承認 URL が必要である。
* `execution:cs10` の接続、`queue.py enqueue` / `cancel` / `pause` / `resume`、reservation の置換、dispatcher / tmux、job の開始・停止は別々の外部操作であり、それぞれに当該操作のユーザー明示承認が必要である。queue reservation は fixed-commit execution contract とし、enqueue 後の branch 移動から置換や cancel を推測しない。
* Phase 2 の `model_profile.implementation_status: pending` には `model-development-evaluation` skill と `development_evaluation` contract を用いる。Issue 固有 evaluator は追加せず、長時間 stage の前に short-screen QC / replay を行う。
* 関連作業は GitHub-native relationship で管理する。bounded child task は parent の sub-issue にし、実際の完了依存だけに `blocking` / `blockedBy` を付ける。
* 新規 Issue では `Roadmap category (Milestone)` を必須とし、`roadmap:triage` / `roadmap:needs-review` の Issue は metadata 修正まで実装しない。
* source Issue を PR から参照する。merge で Issue が完了する場合だけ `Closes #<issue>` を使い、それ以外は残作業を示す。完了には local `review_result.json` の `status: PASS` が必要である。
* file-changing Issue の初回完了報告は、local PASS、commit、push、source Issue を参照する PR 作成後に限る。進捗 commentary は可とするが、PR 作成前に実装完了・成果物・PR 候補を報告しない。PR 不要の明示指定または PR 作成失敗時だけ例外とする。
* merge には required checks と `codex-review-gate` の pass が必要である。物理解釈、dataset 採択、Phase 境界、output contract、ML policy の変更はユーザー明示承認なしに merge しない。merge 後の Issue 引継ぎ手順は `docs/codex/codex_workflow.md` に従う。

## 配置、文書、再現性

* `scripts/` は user-facing CLI / orchestration、`src/` は再利用可能な実装、`conf/` は再現可能な runtime 設定、`schemas/` は machine-readable contract、`docs/codex/` は Codex 運用文書、`tools/codex/` は workflow 補助を担う。
* Phase 文書の正本は `docs/codex/phase_document_policy.md` とする。`phaseX_current.md` は現在地、`phaseX_tasks.md` は判断記録である。統合・移行・削除には `phase-document-maintenance` skill を使い、参照切れを確認する。
* output は JST の `outputs/YYYY-MM-DD/HHMMSS/` に保存し、該当 run では `run.log` と `manifest.json` を残す。Phase 2 は `step_summary.csv` を用い、`step_summary_full.csv` を再導入しない。
* cs10 simulation、archive analysis、render 後は、ユーザー確認用 artifact・manifest・summary・必要な reanalysis archive をローカルへ同期し、件数、SHA-256、QC を確認する。operational log と credential は同期しない。

## 検証と報告

* 最小の targeted test から始める。docs-only、planning-only、workflow-only 変更では full pytest を既定で要求しない。長時間 simulation、sweep、training、render はユーザー明示依頼なしに実行しない。
* Cloud review はユーザー承認後の merge-ready final candidate に対して一度だけ依頼する。actionable review thread は merge 前に解決する。
* 最終報告には、要約、変更ファイル、実行・未実行の check、user review、`review_result.json`、文書と ADR、commit、push、PR、残作業を記載する。
