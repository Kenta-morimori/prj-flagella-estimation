# Issue #255：frame・motor修正候補の検証契約

## 今回の範囲と条件

ローカルの実装・数値検証・実行準備まで。simulation、cs10接続、予約、dispatcher/job開始は実施しない。既存queue #21を操作しない。物理モデルの採択・canonical変更は行わない。

ユーザー指定を更新し、nf=3～5の各1条件、両修正を適用した候補Cのみ、計3条件を用意する。

| 条件 | 付着slot | frame修正 | motor修正 |
|---|---|---|---|
| nf03__slots024 | 0,2,4 | energy_gradient | minimum_norm |
| nf04__slots0134 | 0,1,3,4 | energy_gradient | minimum_norm |
| nf05__slots01234 | 0,1,2,3,4 | energy_gradient | minimum_norm |

A（frameのみ）・B（motorのみ）は設定として選択可能だが、今回のcampaignには含めない。Cと#20を比較してもA/Bそれぞれの寄与は識別できない。単独候補の追加は結果レビュー後の別判断とする。

初期hookは既存の中立・菌体長軸直交配置を維持する。T=2.5e-20 N m/flagellum、dt_star=1e-4、reference τ=0.04 s、局所one-ring全vector反作用、Brownian/switching OFF、seed=0、菌体–べん毛排除OFF、べん毛同士の反発ONを維持する。hook剛性・初期配置の変更と時間発展中の90°拘束は加えない。

## 修正の定義

### frameポテンシャル

`motor.attach_frame_reaction: energy_gradient`を明示選択する。既定値`legacy`は従来の力を維持する。

E = Σ k/(2 l₀²) |(r_right−r_left)−R(r_body) h_local|²。

付着frame Rを定義する長軸（両端layerの重心差）、付着layer重心、付着ビーズのradial方向、正規化と外積の位置依存をすべて解析的に微分する。べん毛・付着ビーズへの従来の直接力に、frameを構成する菌体ビーズへの力を追加する。有限差分は検証にだけ使用する。縮退frameは明示エラーで、固定world軸へfallbackしない。候補はvector tangent modeのみ対応する。

### motor駆動

`motor.segment_torque_correction: minimum_norm`を明示選択する。既定値`none`は従来の駆動力を維持する。

既存segment coupleと実際のdiffusive weightsからF₀を計算し、べん毛ビーズ上でmin ||ΔF||₂、ΣΔF=0、Σ(r−r_root)×ΔF=T a−τ(F₀)を解く。aは既存のPCA駆動軸。補正後の実トルクをone-ringへ逆符号で返す。frameの反作用はframeポテンシャルの微分から求め、motor反作用と混同しない。

この補正は余分な横方向トルクを抑える候補であり、segmentの小さい横方向armによる大きな力を必ず解消するものではない。補正量、補正前後の横方向トルク、補正後駆動力の増大を診断する。solverが解けない場合は失敗し、別方式へfallbackしない。

## 数値検証・診断出力

- frame：有限差分勾配、総力・総トルク、任意の剛体回転・並進、初期ゼロ力、縮退時エラー。
- motor：変形状態・軸に近いsegment、最小ノルム性、軸方向T・横方向ゼロ、総力と反作用、縮退時エラー。
- 既存hex13/project3初期形状と既存方式の非影響、3条件のpreflight、1τ/1s間の物理条件一致をtargeted testで確認する。
- 同一pre-step位置からspring/bend/torsion/hook/frame/repulsion/motor/totalの総力・総トルクを記録する。原点はその時刻の全ビーズ重心。hook欄はframeを除いた成分。
- `force_evaluation_t_s`で力の評価時刻を明示する。既存の形状欄はpost-stepであり、同じ行のt_sと区別する。
- 補正量・トルク誤差、support最小/最大数、solver成功、fallbackをdiagnostic samplesと全step online summaryへ保存し、極値の評価時刻も保存する。失敗はcondition例外として記録する。
- `development_evaluation`のfinite/body/hook length/flag/motor閾値は維持する。hook角度と新規診断欄は判定閾値を追加しない。

## 段階実行と所要時間

campaignは`conf/phase2_multi_run/2010_hex_project_frame_drive_1tau_issue255.yaml`と`2010_hex_project_frame_drive_1s_issue255.yaml`。parallel-jobは`conf/phase2_parallel/issue255_motor_reaction/hex_frame_drive_1tau_job.yaml`と`hex_frame_drive_1s_job.yaml`。

1τは10,000 step、1sは25τ指定で250,000 step。両jobは3条件別output、全条件geometry preflight、`cs10_qualified`・実効3 workers・数値ライブラリthread=1。3独立条件のため8 workersを同時稼働させる追加taskはない。

1τの初期preflight・共通QC・replayをレビューしてから1s対象を確定する。1τでの束化だけでは除外しない。1s用YAMLは候補の準備であり、1τ結果による条件絞り込みや開始の許可を代替しない。新しいcs10操作は固定commitと対象を示して実行許可を確認する。

旧#20の同一3条件の実測を単純比例すると、1τは約16.0/20.9/26.8分、1sは約6.68/8.70/11.18時間。3並列jobの参考wall timeは最長条件相当。ただし旧実行のworker競合、修正solver、追加診断費用が異なるため予約時刻の確約には使わない。Mac見積りは30分超としてheavy simulationを行わない。1s予測は修正後1τの実測で更新する。

## 実行後の比較契約

必要artifactをローカルへ同期し、件数・SHA-256を照合し、run_summaryを先に読む。共通evaluatorと固定cameraの3D/2D replayを用いる。#20は同一条件の0～0.04 s、0～1.0 sに時間範囲を揃え、source manifest・commit・設定一致を確認して対照にする。旧#257の排除OFF結果は初期hook・反作用方式が異なる診断用参考であり、因果比較にしない。

束の広がり・べん毛間軸角・束軸と菌体軸の角度・遊泳方向との角度・速度変動・軌跡の曲がりを、定義とsampling intervalを明記して時系列で比較する。束化は軸角のみで判断せず広がりと併記する。低速で遊泳方向が定義できない時刻は欠測として扱い、角度ゼロに置き換えない。共通QCと遊泳診断を分けて報告し、新しい挙動FAIL閾値は設けない。

Cで改善してもframe/motor単独の効果は未識別である。改善しない場合も横方向反作用だけを原因と断定しない。単独候補、hook compliance、初期geometry、排除体積の追加検証はユーザー判断後に別段階で行う。2sは今回の実行対象に含めない。
