# #263 接続/菌体接触フェーズの実行準備

2026-11-13の新たな自律検証・予約許可に従う。G1未成立の待ち時間に現JのBF ON/OFFを先行し、2条件shortQC/replay→同固定sourceの1秒比較とする。

新規実装はdefault-offの実segment接触テレメトリ。productionの同じ最近点/反発passを観測し、BF/interFF/selfFFの最短距離・作用力・pre-step初接触時刻/segment witness・cutoff内時間を記録する。力/QCは不変。共通summary contractへ単位・欠損・時刻を記載。旧proxy列は残す。

新8targetedPASS、targetgroup60PASS、full1146PASS、e32比30array bitwise一致、Macのshort/main dryrunと共通初期previewを確認。初期dryrunはsingle-key axisのID正規化差で拒否され、linked-key mappingへ修正してYAMLとplanの2ID一致を確認した。simulation予約の失敗ではなく、予約前の設定検査で解消した。

SSH user Ktakemoriと接続成功、queueempty/unpausedを確認。nativeVPNservice表示だけで接続失敗を断定せず、VPN設定は触らない。dispatcherprocessは不在、旧queue tmuxはidle shellで、今回必要な起動は許可済み。旧予約変更・新物理採択・mergeなし。
