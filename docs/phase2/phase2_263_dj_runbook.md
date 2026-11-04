# Issue #263：D/J 2×2の実装と実行手順

2026-11-04追記：ユーザーはR/Dの非重複修正と4条件維持、short QC/replay後の3.5秒予約を承認した。現行の修正/実行手順は[非重複再実行runbook](phase2_263_dj_clearance_runbook.md)を参照。以下の旧結果・判断は当時の固定条件として保持する。

2026-11-03。ユーザーはD/J方針を了承し、cs10実行予約までの作業を許可した。別スレッドでPR #264を中心とする2010 hex modelのフローマップ作成も依頼した。本書は**予約前の実装・準備記録**。追記：予約#40は完走・strict QC4/4 PASSだが、R/Dの初期root対body-edge重複が判明したため3.5秒は未予約。[結果と停止根拠](phase2_263_dj_1tau_result.md)を現状の正本とする。実際のfixed commit・予約ID・開始状況はIssue/PRの要約コメントとローカル予約manifestを正本とする。

## 実装と4条件

| ID末尾 | motor | 初期geometry/自然形状 |
| --- | --- | --- |
| `R_control` | 既存segment couples・拡散重み | 従来compact |
| `D_vector` | `effective_vector_torque` | Rと座標完全一致 |
| `J_transition` | Rと同じ | `axis_transition_rear`、近位3 bend/torsion参照/鎖を事前固定 |
| `DJ_combined` | Dと同じ | Jと座標/自然形状完全一致 |

IDの共通prefixは`nf03__slots024__`。nf3、slots024、FF/BF ON、T/reference T各2.5e−20 N m/本、dt4 µs、粘性1 mPa s、Brownian/切替OFF。解析torsion、指数型反発・修正距離・材料係数を維持する。既定profileは変更していない。

Dは現行と同じ全鎖SVD軸に対してゼロ合力・純vector T eの力をゼロから生成し、one-ringへ逆vectorを返す。拡散状態はDでは時間発展も力生成への使用もしない。縮退solverにfallbackしない。Jはroot半径0/.035/.040/通常半径、hookを外向きから後方10°にし、結合長を保持する。近位自然角はimmutableなmodel参照値としてengine・共通step-summary QC・preflightで共有し、時間発展中の形状へ追従させない。自然角の変更をQC閾値の緩和で代用しない。

`initial_helix_axis_from_rear_deg:0`について、Jは構成したscrew軸が後方を向く。全鎖SVD軸は別診断で、約1.5°の差を記録する。既存compactの「SVD軸誤差ほぼ0」という幾何前処理条件を新しい自然形状へ誤適用しない。物理QC閾値は不変。

## 初期条件の成立と比較の限界

- 現行/Dの束半径指標0.813333 µm、J/DJは0.761913 µm。開始条件Cの1.040678 µmより両群とも小さく、従来compactを対照にする近束screenとする。**初期から完成した密束ではない**。密束維持の主張は、実際に束の受入指標を満たす時刻以降だけに限定する。分離配置からの自発束形成を成功条件に追加しない。
- hook0.25 µm、各bond0.58 µm、bead数を保持。初期弾性力ノルムは現行/D 2.04e−26 N、J/DJ 2.64e−26 N。current/Dの初期菌体反発ノルムは1.039 pN、J/DJは0。これもJの配置変更の影響として示す。
- Jの全位相・未変形body-clearance下界は約0.211 µm（比較距離0.20）。body中心凸包を包含する円筒への距離を位相について最小化し、線分上の有限samplingをLipschitz補正する。動的非貫通は別の診断が必要。
- 遠位の通常らせん形状・物性は共通だが、遠位の絶対座標・位相表現は両geometry群で完全一致しない。Jは根元の自然形状と初期配置を含む構造対照で、「接続点だけ」の純粋な因果効果とはしない。Dの対比較は各geometry内で座標一致する。
- productionの受動Hessian/RPY検査は、現行dt*lambda=0.102、J=0.850。差分幅1/0.5 pmの両方で一致。motor/接触を含む非線形精度・動的QC・遊泳成立を保証する結果ではない。

## 予約・実行順

| Stage | campaign | parallel job | 時間・step/条件 |
| --- | --- | --- | --- |
| short | `2010_hex_nf3_dj_1tau_issue263.yaml` | `nf3_dj_1tau_job.yaml` | 1τ、0.04秒、10000 |
| main | `2010_hex_nf3_dj_3p5s_issue263.yaml` | `nf3_dj_3p5s_job.yaml` | 87.5τ、3.5秒、875000 |

campaignは`conf/phase2_multi_run/`、jobは`conf/phase2_parallel/issue263_initial_bundle/`。どちらも独立4条件、`max_workers:auto`/`cs10_qualified`（実効4 workers）、数値thread1、条件別output、全条件geometry preflight。simulationなしdry-runをMac/cs10で確認してから予約する。

既存dispatcherを利用し、既存予約の変更・取消・置換はしない。dispatcherが不在の場合のみ、今回承認されたjob開始のために起動する。予約はfixed commitの契約であり、その後の記録commitで置換しない。

長期を先に無条件enqueueしない。短期完走後はportableを同期して集合/件数/全SHAを確認し、run_summaryを先読、共通strict QC・固定camera 3D/2D replayをレビューする。明白な崩壊・貫通/実装不具合・QC FAILならlongへ進まない。単に1τで十分な周期数がないことは長期を止める理由としない。短期の成立確認後、合意済み3.5秒4条件を予約する。

旧nf3 ON実測からshort約40分（共有負荷等を含め約40〜60分）、main最遅約55時間が予約前目安。D/Jの実測checkpoint速度で更新する。Macで長時間simulationは実行しない。

## 検証と採否

Dの全vector/仕事/単位/回転並進整合/縮退、root reference共有、初期幾何、長短時間の条件一致、32内部stepの限定積分をtargeted testで確認した。full pytest1108 PASS。共通QCはfinite/body/hook length/flag/motorを維持し、hook対body角は診断専用のまま。

長期の機序窓、緩和後flag/body周期数、束/回転/前進、摂動・数値感度の受入案は[修正方針](phase2_263_model_repair_direction.md)を使用する。既存の原因監査も保持する。短期予約は科学的採択ではない。新しいhook構成則・駆動軸・硬い接触則を追加採用する場合は根拠と比較範囲を示す。

準備正本：`outputs/2026-11-03/183200/issue263_dj_preparation/`。geometry_preflight/初期座標NPZ/preview、production_passive_stability、targeted/full test log、長短dry-run、manifest/run.log/review_resultを保存する。予約の証拠は同日JSTのreservation outputへ分離する。追加位置復元処理/scaleの削除は従前のmerge前残件で、今回の全条件では無効（scale1）。
