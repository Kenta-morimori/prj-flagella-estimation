# Issue #263：反発OFFの束回転鈍化監査

- 主checkoutは別タスクのbranchなので、PR #264の既存worktreeで診断文書のみ編集。
- #35の共通strictQC PASS・全件照合済みローカルarchiveを再利用。run_summaryを先読し、全13 OFF条件の全べん毛spinを走査した。
- 5条件各177状態をdt=0で力評価し、各力とmobilityの線形寄与を算出。微小仮想変位で角速度演算子を確認。物理時間の更新や動的な力OFF試験は行っていない。
- nf6では近位segment接触、全鎖の同時減速、接触解放後の回復を確認した。
- 線分間距離の単位依存を静的に再現。cutoffの不連続も式から確認したが、軌道への寄与の分離は未完了として残した。
- 結果と3条件の比較案を日本語文書に整理。相互作用修正は未着手・未予約でユーザー判断待ち。
- 新たなQC evaluator、torsion方式比較、モデル採択、mergeなし。監視automationはPAUSEDを維持。
- 専用の一時venvを作り、他タスクの共有editable environmentを変更せず既存hookを実行する。
