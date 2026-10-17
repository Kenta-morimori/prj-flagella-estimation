# Issue #255：局所反作用2s結果（逐次追記）

この文書はqueue #20（2010 hex）とqueue #21（2010 project）を別々の完了単位として記録する。2026-10-08時点では**#20のみ評価済み**。#21の予約・実行・成果物には触れていない。両jobの実行commitは`e94c7cd4f0e31f99169b28ff5025e58e3b02b7cc`。物理モデルやcanonical profileの採択は行わない。

## #20：2010 hex 13条件

cs10のqueue #20は2026-10-06 17:09:50 JSTに開始し、2026-10-08 10:09:29 JSTに完了した。13/13 conditionが各500,000 step、2.0 sまで完走し、parallel aggregationも完了した。cs10のcampaignからportable artifact **122ファイル**をローカルの`outputs/2026-10-08/171556/issue255_hex_local_2s/`へ同期し、リモートとローカルのファイル集合・SHA-256を全件照合した。operational logは同期していない。ハッシュ一覧は同じ出力ディレクトリの`sync_verification.json`に保存した。

`development_evaluation`のlong-duration判定は**13/13 PASS**。全stepのfinite・body・hook長・flag・motorゲートにFAILはなく、`run_summary.json`の行数・step連続性も全条件で一致した。最大motor力残差比は`3.97e-15`、トルク残差比は`4.36e-15`で、契約上限`1e-8`、`0.02`を十分下回る。motorのdegenerate axis / split rank deficient countは全条件で0。hook長相対誤差の全条件最大は`0.1002`、flag bond相対誤差最大は`0.0954`、body spring stretch最大は`0.0912`だった。全条件の反作用設定は`attach_one_ring`で、保存されたtopologyでは各べん毛のsupportは付着ビーズを含む5ビーズ。局所solverが失敗すればこの実装は例外でconditionを失敗させるため、13条件の完走とmotor残差は局所solver成功を支持する。**per-stepのsolver/support/fallback telemetryはartifactにない**ため、観測値としてfallbackゼロとは記さない。設定上のfallback modelは`none`。

### 形状診断

初期geometry preflightでは全13形状の後方軸角0°、hook平衡角誤差最大`1.5e-14°`、hook長誤差最大`1.1e-22 m`、初期hook曲げ力最大`9.1e-28 N`、attach→firstと菌体長軸の直交誤差0°だった。ビーズ非重複・外向き配置もpreflightを通過した。

時間発展中のattach→first対菌体長軸角は**診断専用**であり、90°からの偏差をFAILには使っていない。全step最大偏差は`65.93°`（`nf03__slots012`）、保存された約1 ms間隔の2,002状態による角度時系列は`diagnostics/hex_2s_hook_axis_timeseries.csv`と図に保存した。初期値はほぼ90°だが、2s中に大きく変わる条件がある。

|条件|本数|QC|hook長相対誤差max|flag bond相対誤差max|後方軸角終端|hook長軸偏差max|
|---|---:|---|---:|---:|---:|---:|
|`nf01__slots0`|1|PASS|0.092|0.060|69.2°|33.1°|
|`nf02__slots01`|2|PASS|0.100|0.074|73.2°|37.4°|
|`nf02__slots02`|2|PASS|0.095|0.064|61.8°|37.8°|
|`nf02__slots03`|2|PASS|0.094|0.071|48.1°|27.1°|
|`nf03__slots012`|3|PASS|0.100|0.089|103.7°|65.9°|
|`nf03__slots013`|3|PASS|0.094|0.088|73.6°|41.5°|
|`nf03__slots014`|3|PASS|0.096|0.086|117.5°|42.2°|
|`nf03__slots024`|3|PASS|0.092|0.071|57.5°|39.2°|
|`nf04__slots0123`|4|PASS|0.100|0.085|120.6°|50.1°|
|`nf04__slots0124`|4|PASS|0.097|0.092|129.6°|40.7°|
|`nf04__slots0134`|4|PASS|0.094|0.077|49.2°|41.2°|
|`nf05__slots01234`|5|PASS|0.099|0.088|79.6°|54.7°|
|`nf06__slots012345`|6|PASS|0.098|0.095|100.1°|51.7°|

後方軸角の終端値が90°を超える条件は5/13あり、初期の後方整列が2s後にも保たれるとは言えない。らせんpitch相対誤差の全step最大は`1.28e6`で、pitch fitの不安定な外れ値を含む。終端でも最大`28.85`（`nf02__slots01`）。pitchは現行QCの判定項目外であり、この値だけで物理的な巻きの崩壊を断定しない。固定cameraの3D終端画像には大きく曲がるべん毛が見えるため、形状解釈は動画・保存座標と合わせた追加レビューが必要。

### ローカル成果物

- 共通QCと12本の固定camera動画（n=1～6ごとに3D/2D各1本、各41フレーム）：`outputs/2026-10-08/171556/evaluation_hex/`
- 初期2D配置、全条件のQC表・heatmap：同ディレクトリの`initial_geometry/`、`summary.csv`、`heatmaps/`
- 軸方向の初期/終端2D図、角度時系列図、13条件の定量CSV：`outputs/2026-10-08/171556/diagnostics/`
- 元archive・run summary・provenance、同期ハッシュ：`outputs/2026-10-08/171556/issue255_hex_local_2s/`、`outputs/2026-10-08/171556/sync_verification.json`

共通2D動画は主に菌体輪郭を表示し、べん毛形状の判定には軸方向2D図と3D動画を併用する。図は固定座標系で比較し、菌体の移動・回転も残している。

## #21：2010 project n=1～3（未完・未評価）

job完了後、#20とは独立したローカルcampaignディレクトリへ同期・SHA-256照合し、同じ`development_evaluation`のproject用configでQCとreplayを作る。この節に3条件の実測、定量表、形状所見を追記してから16条件の総合判断を行う。現時点で#21の結果を推定して埋めない。
