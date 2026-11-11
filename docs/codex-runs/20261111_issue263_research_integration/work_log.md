# Issue #263：research結果の統合と修正候補の静的検算

2026-11-11。別researchチャット完了後、ユーザーの依頼に基づき親側で結果を統合した。Issueのexecution:cs10 / cs10_user_run・milestone2とPR264を確認。primary checkoutは別タスクのため変更せず、既存PR専用worktreeを使った。

researchのreport/sourcesとmanifestをSHA照合し、新情報Bianchi2023本文/補足の式と近似を親側でも確認。既存#42のrun_summaryを先読し、portable37件の集合・件数・SHAを再確認。src/conf/testsは実行59890c0と同一。

局所接続を一度に新設する案から、現D+Jと同じ座標・支持点で相対位相の仕事を導くM_phaseへ具体化。位置差分prototypeで保存13状態×3鎖の力・torque・仮想仕事、入力/出力spin・純tilt、客観性、phase gauge、反平行拒否を検算した。位相testの219点と、独立剛体spin全域の非重複証明は別に扱った。後者は現Jの固定軸回転族だけに成立し、密束Gは未成立のまま。

bodyを含む全RPYとhook長の瞬間拘束を使う4剛体群referenceを作り、仕事と正の散逸を確認。初期はMの首振りrateが小さいが、後半保存状態では増える鎖もあった。反実仮想の瞬間応答を実軌道の改善予測にせず、反証も記録した。

次の案はDJ/M/K/MKの2×2。出力軸T、dt、粘性、J座標・自然形状、FF/BF等を固定する。全vectorの大きさ・仕事率はMで変わることを明示。production解析Jacobian、共通診断/targeted test、.04秒QC/replay、約1秒機序比較、近束維持の入口を具体化し、実装/採択/予約はしていない。静的prototypeはignored outputsへ保存し、Issue専用のacceptance evaluatorを追加していない。

新規予約保留は維持。既存メール通知設定・他chat・cs10 queueは変更しない。ADR0029はProposed、追加位置復元code/scale削除はmerge前残件。重要な結論をIssue既存要約、PR本文/既存要約へ反映し、詳細は設計正本へリンクする。
