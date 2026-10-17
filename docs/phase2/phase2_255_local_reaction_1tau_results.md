# Issue #255：局所・全ベクトル反作用の1τ結果

2026-10-05、固定commit `2e8f91a6913e0ef22de58337cfa0aaaf63c9175a`をcs10で実行した。queue #18は2010 hexの13条件、#19は2010 projectのn=1..3の3条件で、両jobとも全condition成功した。project n=4..6は対象外である。1τは10,000 step、0.04 s。hexは42分41秒、projectは9分41秒、連続した両jobは約52分22秒を要した。

成果物はローカルの`outputs/2026-10-05/210315/`にある。hex 122ファイル、project 32ファイルをremoteとローカルで件数・SHA-256照合した。子run manifest 16件も個別にSHA-256照合し、局所supportは全件で付着点とone-ringの5ビーズ、fallbackは全件falseだった。成果物全体のローカルSHA-256一覧は`diagnostics/manifest.json`、図・動画への入口は`diagnostics/review.html`である。operational logは同期していない。

## 共通short screen

`development_evaluation`の全step判定。PASS数はfinite、body、hook length、flag shape、motorのすべてを通過した条件数。既存の軸方向・全body結果は固定commit `6592e140e225ea1ce981700d6854bcc99ab707e4`の1τ成果物を同一形状の対照に用いた。

|モデル|菌体側反作用|全ゲートPASS|最大motor torque残差|最大hook長相対変化|最大flag bond相対変化|最大\|hook角−90°\||
|---|---|---:|---:|---:|---:|---:|
|hex|軸方向|4/13|0.2963|0.03586|0.02929|11.03°|
|hex|全body・全vector|13/13|3.62×10⁻¹⁵|0.03574|0.03241|11.69°|
|hex|局所・全vector|13/13|2.80×10⁻¹⁵|0.03936|0.03674|12.22°|
|project|軸方向|1/3|0.1960|0.04206|0.05165|0.72°|
|project|全body・全vector|3/3|3.37×10⁻¹⁵|0.04209|0.05336|0.94°|
|project|局所・全vector|3/3|1.25×10⁻¹⁵|0.04459|0.05634|1.18°|

局所方式のfinite/body/hook length/flag shape/motorは個別にも各16/16 PASS。最大motor force残差も数値精度の範囲だった。局所solver失敗・全body fallbackは0件。16形状ごとの値は`diagnostics/three_arm_comparison.csv`に記録した。反作用の釣り合いは改善したが、この結果だけで局所方式の物理的採択は確定しない。

## 角度・形状・replay

付着body bead→第1flagellum beadと菌体長軸の角度はt=0で90°。専用角度ポテンシャルはなく、角度はQC閾値ではない。上表の最大偏差はrun summaryから**全10,000 step**を対象に算出した。局所方式の各べん毛の時系列は保存archiveの41時点を`diagnostics/local_hook_axis_timeseries.csv`とモデル別plotに記録した。hexは最大12.22°で、終端まで増加傾向が続く。projectは最大1.18°である。

固定cameraの3D/2D replayは9組・計18動画、各41フレームを確認した。3Dではらせんの大きな破綻・崩壊は認めなかった。汎用2D replayは菌体が主に写り、べん毛形状の判定に向かない。このため保存した3D座標から、初期と1τ終端を同じ菌体長軸方向に投影した`diagnostics/hex_axial_initial_final.png`、`diagnostics/project_axial_initial_final.png`を補助図として作成した。両図でも全形状のらせんが判別でき、大きな崩壊は見られない。局所方式は全body方式に比べ、最大hook長・flag bond変化とhexの角度偏差がやや大きい。

## 2sへの判断材料

局所方式の16形状は1τの数値・形状screenを全件通過し、2s候補になり得る。一方hexの角度偏差が1τ終端で増え続けており、長時間安定性は未確認。局所方式の16条件を同じ8/3 workersで50倍のstep数へ延長する単純比例見積りは約43.6時間で、実際の所要時間は未検証である。2sの対象と開始はユーザーレビュー後に決める。2sは開始していない。ADR 0024はProposedのままで、canonical modelも変更していない。
