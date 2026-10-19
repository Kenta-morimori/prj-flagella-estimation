# Issue #268 work log

- Source Issue: https://github.com/Kenta-morimori/prj-flagella-estimation/issues/268
- Parent Issue: https://github.com/Kenta-morimori/prj-flagella-estimation/issues/266
- Branch: `codex/issue268-sol-harness`。`gh issue develop` で GitHub Development に連携した。
- Execution target: `no_runtime` / `execution:none`。独立 condition と Mac heavy wall time は N/A。simulation、sweep、training、render、cs10 操作は行っていない。

## 第1段階: モデル移行

- `.codex/config.toml` を `gpt-6.1-sol` / `medium` に変更。TOML parse で値を確認した。
- OpenAI Docs の `gpt-6.1-sol` model page は `medium` をサポートする。アプリ同梱 CLI `0.162.0-alpha.17.2` の短い read-only 起動でモデル・effort・応答を確認した。
- ユーザー判断で旧モデルとの事前比較を省略。成功率、所要時間、credit の改善は主張しない。旧 CLI `0.148.0` はこのモデルをサポートせず、評価にはアプリ同梱 CLI を使った。
- 参照検索では旧モデル名は初期設定の判断履歴と今回の例外説明にのみ残した。

## 第2段階: 指示評価

- `instruction_audit.md` に分類、正本、変更候補と根拠を記録。`evaluation_cases.md` と `eval_prompts/` に7ケースの期待結果・禁止操作・再実行可能な prompt を保存した。
- ユーザーの外部送信承認後、現行版と Phase 文書整理 skill ルート短縮候補を同じモデル・effort・read-only sandbox で比較。候補全文を `candidate_phase_document_maintenance_SKILL.md` に保存し、`evaluation_results.md` に14回の結果を記録した。
- 両版とも7/7ケースで期待結果と安全境界を満たしたが、候補は Phase 文書ケースの入力 token・command 数を減らさなかった。短縮と skill の移動・symlink は見送り、AGENTS・skill・policy test を変更しなかった。
- cs10 操作別承認、Phase 2 QC、local PASS と CI / review gate の分離を保持した。料金・credit は取得できず、所要時間はファイル時刻の代理指標とした。

## 月次運用と検証

- `docs/codex/monthly_harness_review.md` を #252 索引、月次親、採択事項の子 Issue、親の完了条件へ更新。今回の事前比較省略を例外として記録し、通常の比較手順を維持した。ADR 0029 に判断を記録。
- 対象 policy test 6件、TOML 構文・設定値、`git diff --check`、旧規則・参照検索を確認。docs/config/workflow の変更であり、simulation と full pytest は実行していない。
- Local review: 物理モデル・dataset・output contract の変更なし。参照、正本、Issue #266 の判断との整合を確認した。
