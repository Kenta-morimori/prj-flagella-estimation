# Issue #268: 同一モデル指示評価結果

2026-10-19 にアプリ同梱 Codex CLI `0.162.0-alpha.17.2` を使い、`gpt-6.1-sol` / `medium`、read-only sandbox で `evaluation_cases.md` の7ケースを現行版と候補版に1回ずつ実行した。候補は Phase 文書整理 skill のルート文書だけを231行から25行へ縮め、詳細は既存 checklist と policy に案内した。使用した25行の全文は `candidate_phase_document_maintenance_SKILL.md` に保存した。現行版の元ファイル `tools/codex/skills/phase-document-maintenance/SKILL.md` の SHA-256 は `b02834b5439b68f500b90e47acbb58974aae7c9b159d842069704f938b3ffd73`、候補の SHA-256 は `6a142eb21c9c1762a50b9f2bc03218dae7281a7dc48958be7fa81de3a942f990`。比較対象の repository コピーでは、それ以外の内容と prompt を同一にした。ユーザーはこの評価で repository 内容を Codex/OpenAI へ送信することを明示承認した。

| 指標 | 現行版 | 候補版 |
| --- | ---: | ---: |
| 期待結果を満たしたケース | 7/7 | 7/7 |
| 禁止操作・安全境界違反 | 0 | 0 |
| Phase skill 選択（適用／非適用） | 1/1、非適用6/6 | 1/1、非適用6/6 |
| 入力 token（cache 込み、CLI報告値の合計） | 540,168 | 616,989 |
| cache 済み入力 token（内数） | 415,872 | 489,344 |
| 出力 token | 4,407 | 4,598 |
| command 実行数 | 15 | 19 |
| JSONL ファイル作成から最終更新までの合計秒数（所要時間の代理指標） | 143.2 | 137.5 |
| Phase 文書ケースの入力 token／command 数 | 92,571／3 | 123,872／4 |

両版とも Phase 文書ケースで skill、正本 policy、checklist を読み、判断保持・live contract・削除 gate・参照確認を説明した。非適用ケースで Phase skill を誤選択せず、Issue 化、pending model、cs10、PR 引継ぎの境界も保持した。全 command は読み取りで、repository・GitHub・cs10 を変更していない。code ケースはどちらも実装場所と対象 test を示したが、候補版では調査 command が3件から6件へ増えた。

**判断:** 候補の短縮は採用しない。対象 Phase 文書ケースを含めて token と command 数の改善が確認できず、秒数の差は単回実行かつファイル時刻に基づく代理値なので優位性の根拠にならない。現行の `AGENTS.md` による明示導線で適用ケースは skill を見つけ、非適用ケースでは発火しなかったため、`.agents/skills` への移動・symlink も今回は見送る。指示の意味を変える編集や文字列重複を強制する policy test の変更は行わない。

token は CLI の報告値で、料金・credit は取得していない。単回の小規模ケースなので、成功率や本番 task の所要時間改善は推定しない。次回月次親で移行後の実績を確認する。旧 CLI `0.148.0` は `gpt-6.1-sol` を未知モデルとして扱い、ChatGPT アカウントで非対応の API 応答となった。アプリ同梱 CLI の短い通常 task は正常に起動した。
