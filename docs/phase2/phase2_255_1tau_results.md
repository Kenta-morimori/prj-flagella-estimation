# Issue #255：1τ short screen 結果

2026-10-05。対象は hex 13配置×2方式、project n=1..3×2方式、計32条件。
cs10予約 #16/#17 はともに succeeded。simulation固定commitは
`6592e140e225ea1ce981700d6854bcc99ab707e4`。各条件10,000 step、1τ=0.04 s。
2sは未実行であり、対象選定と開始はユーザーレビュー後に判断する。

## 成果物と検証

ローカル成果物基点は `outputs/2026-10-05/142734/`。
`issue255_1tau_hex/` の239ファイルと `issue255_1tau_project/` の59ファイルを
既存sync helperで同期し、remote/localの件数とSHA-256を照合した。
`run_summary.json` を先に確認し、共通development evaluatorで評価した。
`diagnostics/manifest.json` に成果物ハッシュ、`diagnostics/review.html` に映像一覧を保存した。
`evaluation_hex/`、`evaluation_project/` にQC、初期2D配置、固定camera replayを保存した。
20本のMP4をdecodeして各41frameを確認し、全32条件の中間・終端3D画像を目視確認した。
確認した画像に明らかな菌体崩壊やべん毛形状の破綻はない。これは1τ内の確認に限る。

## QC

|判定|hex|project|合計|
|---|---:|---:|---:|
|finite/body/hook length/flag（各項目）|26/26|6/6|32/32|
|motor：軸方向反作用|4/13|1/3|5/16|
|motor：全vector反作用|13/13|3/3|16/16|

閾値と物理モデルは変更していない。motor FAILを含む11条件は診断専用。
全vectorのmotor torque residual最大は3.63e-15未満、軸方向の最大はhex 0.29632、
project 0.19601（既存閾値0.02）。条件別値は `diagnostics/qc_breakdown.csv`。

## hook角度

attach→firstと菌体長軸の最大 |角度−90°| はhex **11.69059°**、project **0.94132°**。
最大値は全10,000 stepのonline集計。各べん毛の時系列CSVとPNGは
`diagnostics/hook_axis_timeseries.csv`、`hex_hook_axis_timeseries.png`、
`project_hook_axis_timeseries.png`。時系列は41保存時刻（1 ms間隔）であり、全stepではない。
hexは1τ終端でもずれが増加しており、角度が平衡化したとは判断できない。
角度は診断値とし、専用ポテンシャルや新しいQC閾値を導入していない。

## 2s候補

形状ゲートは全16対が通過。両反作用方式でmotorを含めPASSした5対は次の通り。

- hex：`nf02__slots03`、`nf03__slots024`、`nf04__slots0134`、`nf06__slots012345`
- project：`nf03`

残る11対を比較対象に含める場合、軸方向motor FAILの2s結果は診断専用とする。
`diagnostics/pair_candidates.csv` に候補一覧を保存した。project n=4..6は含めない。
実測wall timeはhex 63.06分（8 workers）、project 9.97分（6 workers）。
全32条件を同じjob構成で2s=50τへ単純比例すると約60.85時間。
これは長時間の性能・安定性を保証する見積りではなく、対象確定後に再見積りする。
canonical model採択、2s開始、mergeは行っていない。

## replay表示修正とレビュー

共通replayが従来の形状判定のみを表示し、motor FAILをPASSと表示する問題を修正した。
development evaluationのstage/statusをreplay入力へ渡し、映像の判定表示へ反映する。
simulationは再実行せず、修正版で全replayを再生成した。
対象テスト81件PASS、ruff check/format、git diff --checkを確認。
変更は解析表示と結果文書であり、新しいADRは不要。初期形状の判断は既存ADR 0023を維持する。
