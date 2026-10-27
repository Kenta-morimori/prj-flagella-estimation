# Issue #263：距離計算修正と再実行準備

- ユーザーは準備を依頼し、予約は許可後と指定。cs10接続・queue・dispatcher・simulation操作は行わない。
- PR #264専用worktreeを継続使用し、主checkoutの別タスクには変更しない。
- 幾何27ケースを先に追加し、修正前9 FAILを確認した。
- 単位に依存する絶対閾値を除去し、正規化・外積による内点解・4境界最小化・縮退処理へ修正した。cutoffや係数は変更しない。
- 小配列の外積allocationをscalar triple productへ置き換え、同じ回帰testを維持しつつ評価費用を抑えた。Macの旧/新median費用比は約1.027。
- nf6修正4us/2usの1tau・25tau設定を用意。旧対照は#35の0〜1秒を再利用し、旧solver切替optionを追加しない。
- 共通preflight・preview・両parallel dry-runを確認。座標は同じbuilderの旧設定と完全一致し、cs10 archiveとの差は丸め程度。初期反発0・弾性力は丸め誤差。
- targeted 55 PASS、既存Issue263初期配置test 8 PASS。共通previewを目視確認した。
- 初回screen約1時間、1秒比較約26時間（2workers並列）の暫定外挿を保存。動的結果・新規cs10 benchmarkではない。
- 日本語runbookとlocal review_resultを追加。モデル採択、torsion追加比較、Cloud review、mergeなし。
