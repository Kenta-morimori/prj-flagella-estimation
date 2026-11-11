# Issue #263：M/Kの解析力と4条件の実行手順

2026-11-11。Refs #263 / PR #264。[承認済み自律ロードマップ](phase2_263_autonomous_roadmap.md)の第1フェーズ。[research統合・力学仕様](phase2_263_connector_research_synthesis.md)に基づくopt-in比較を実装した。科学的採択・既定profile変更は含まない。

**実行状況**：short予約#43を固定`e32c4ad033411c4bbdfb3fcc6583eab10b93e08e`で作成し、11/11 23:02 JSTに4条件running/failed0を確認。[予約レビュー記録](../codex-runs/20261111_issue263_mk_short_reservation/review_result.json)。完走/QC/遊泳成立は未判定で、1秒stageはshortの代理レビュー後。予算は1/5予約使用。証拠は`outputs/2026-11-11/230216/issue263_mk_short_reservation/`。

## 比較と固定値

| 条件ID末尾（prefix `nf03__slots024__`） | motor.force_distribution | potentials.bend.root_response |
| --- | --- | --- |
| DJ_control | effective_vector_torque | cosine |
| M_phase | relative_phase_work | cosine |
| K_internal | effective_vector_torque | normalized_tangent |
| MK_combined | relative_phase_work | normalized_tangent |

全条件で現在Jの座標・固定自然角・局所位相offset0を共有する。nf3/slots024、出力軸T/reference T=2.5e−20Nm/本、dt4µs、粘性1mPa s、FF/BF ON、Brownian/切替OFF、解析torsion、ばね/遠位材質/反発曲線・strict QCを維持する。近束Gや#271配置は取り込まない。現Jは完成した密束ではなく、この4条件は機序診断である。

## Mの解析勾配

`phase_drive.py`はbody局所5 beadとflag全11 beadの固定material templateからproper fitを作る。flag側は現在の向き付きPCA軸へ最小回転で合わせる。相対quaternionの位相χに対し `F=T∇χ` を解析的に計算する。新しい材料twist状態や初期位置への復元energyはない。

raw fitのflag軸a、PCA軸e、body入力軸nを使い、

```text
g = (n+e)/(1+n·e)
k = (a+e)/(1+a·e)
h = (a×e)/(1+a·e)
dχ = k·dθ_fit + (g×e−h)·de − g·dθ_body
```

からproper-fitとPCAの両方へpullbackする。proper fitは、無次元化した `H=YᵀX` のproper polar rotation Rについて `J=tr(HRᵀ)I−HRᵀ` を用いる。couple uに対し `w=J⁻¹u`、各点の勾配は `w×(RX_i)/length_scale`。PCAは主固有vectorと他固有vectorの固有値差を使う解析微分である。変形した鎖で指定coupleを最小ノルム力へ置くだけの実装ではない。

剛体変位ではflag coupleは `Tg`、bodyはその反作用になり、出力軸投影Tは固定される。全vectorの大きさと仕事率は変わる。proper fit非一意、PCA縮退/向き不明、最小回転または入力/出力の反平行は数値guard1e−8で拒否する。clipや別modeへのfallbackはない。内部stepのχ増分が0.25rad以上なら時間分解不能として停止し、dtを自動変更しない。

## Kの保存力

近位の固定自然角を持つ全9 triplet（各鎖の開始0/1/2）だけ、`kb_i=kb*sin(normal_rest_angle)²/sin(root_rest_angle)²` とする。既存cos型energyの勾配を係数配列で評価し、自然角・遠位・bodyは変えない。normal=142°で、倍率は約27.9 / 2.18 / .594。自然角のsin²が1e−8以下では拒否する。全体tiltを戻す力ではなく、近位内部応答の因子である。

## 診断の意味

- `motor_drive_contract_version=1` は現D。旧 `motor_drive_axial_error_ratio_max` は従来通り指定**全vector**との差。名前にaxialがあっても投影誤差へ意味を変更しない。
- version2はM。旧D専用metricはNaN（CSV）/null（online summary）として未適用を明示。物理有限性gateと混同しない。共通strict QCは既存の合力/合torque残差を維持する。
- 両driveの `target_vector_error_ratio_max` と `output_projection_error_ratio_max` を区別する。Mは入力軸投影、入力/出力角、vector gain、fit conditioning、PCA gapと反平行余裕も記録する。phase専用値はversion2で解釈する。
- `phase_work_step_J = Σ T Δχ`、`phase_linear_work_step_J = Σ F_pre·Δq`、累積workと差を診断する。仕事は実step終了時刻、力/torque指標は `force_evaluation_t_s` のpre-stepへ対応させる。仕事診断は力へfeed backしない。
- χは有効接続のspin/precession位相で、材料twistやべん毛自転そのものではない。既存の各鎖材料点spin・菌体材料点自転を使い続ける。有限長hookのFENE/曲げの力と仕事は既存の受動成分として残し、motorへ二重加算しない。

## 実行前に確認した証拠

正本：`outputs/2026-11-11/214352/issue263_mk_preparation/`。

- 独立した位相oracleの任意非剛体変位/全座標差分、客観性・単位共変性・位相gauge、入力/出力投影、1回転の仕事、逆軸/縮退停止、Kのエネルギー勾配、共通QC/テレメトリをtest。
- 既存#42のrun_summary先読とportable37件の集合/件数/全SHAを確認。13保存状態で独立方向の仕事差は最大約1.06e−8（T/µmで正規化）、指定couple誤差は丸め誤差。
- 固定baseline `59890c05460ba3032d50813afaa488ed3440f271` のsrc/confを隔離抽出。D+Jについて初期・64step終端の座標/力、13保存状態のmotor/total力を新実装と比較し、各配列がbit単位で一致した。
- Macのbounded64step費用はDJ .108秒/step、M .110、K .113、MK .122（他check並行中の粗い値）。物理時間は各.000256秒のみ。桁違いの費用増加はないがcs10 ETAの代用にしない。
- Mac dry-runは4条件・実効4workers/thread1・条件別output・全初期active線分clearance PASS。初期プレビューは同一形状1種類を表示し、4条件の同一座標hashを保持する。

全testの詳細と最終commit/入力SHAはreview_resultと予約manifestへ記録する。単体testや64stepだけで1τ QCまたは遊泳成立を宣言しない。

## 固定commitと実行順

| Stage | campaign | parallel job | 各条件 |
| --- | --- | --- | --- |
| short | `conf/phase2_multi_run/2010_hex_nf3_mk_1tau_issue263.yaml` | `conf/phase2_parallel/issue263_initial_bundle/nf3_mk_1tau_job.yaml` | .04秒 / 10000 step |
| mechanism | `conf/phase2_multi_run/2010_hex_nf3_mk_1s_issue263.yaml` | `conf/phase2_parallel/issue263_initial_bundle/nf3_mk_1s_job.yaml` | 約1秒 / 公称250000・runtime250001 step |

実行sourceは本実装の固定commitとし、enqueue前に具体的なSHA・job/campaign SHA・4ID・runtime計画stepを予約manifest/Issue要約へ示す。1秒の浮動小数点ceil差はplanと実終端で照合する。将来branch headが進んでも予約を変更しない。

Macの1秒planは、reference τ=.039999999999999994sと既存ceil実装により250001step・終端公称1.000004sとなる。contractのexpected_stepsも250001と明示し、主比較窓は1.0秒までとする。T/dtや既存step計画の実装を変更して差を隠さない。shortは10000stepのまま。

1. Macのtargeted test、source回帰、preflight・preview・simulationなしparallel dry-runを確認しcommit/push。
2. 許可済みcs10接続でqueueと既存dispatcherをread-only確認。固定commitを取得し、隔離sourceで同じjobのdry-runを行う。既存job/予約の条件を変えない。
3. `queue.py enqueue --branch <fixed commit> --config conf/phase2_parallel/issue263_initial_bundle/nf3_mk_1tau_job.yaml`。FIFO/既存dispatcherを利用し、必要な起動/今回job開始は承認済み。
4. 予約ID/commit/入力SHA/開始正常性・ETAを記録。台帳は`outputs/2026-11-11/212021/issue263_autonomous_roadmap/execution_ledger.json`。enqueueの応答が不明なら同じ予約の状態を照合し、重複enqueueしない。失敗実行も枠に算入する。
5. 完走時はportableのみ同期し集合/件数/全SHAを確認、run_summary先読、共通development_evaluation strict QC/heatmap/固定camera3D/2D replay・frame/fps/coverage/SHAと視認、QC-PASSの束/回転/菌体材料点自転/前進/根元診断。phaseχを材料点spinと混同しない。
6. 代理レビューでshort PASSかつ明白な異常なしなら、同じfixed commitの1秒stageを予約する。周期不足だけでshortから止めない。最大2回のshort→1秒、入口を満たす候補＋DJの3.5秒1回までという台帳上限を守る。

今回の比較窓は`.04–.20 / .20–.65 / .65–1秒`。指標・進行/棄却基準は設計正本とロードマップに従う。未達候補の時間延長、QC閾値緩和、T/dt/粘性の調整、無根拠な再実行はしない。新しい構造変更や3.5秒超の実行へ自動で拡大しない。科学的採択/merge/canonical freezeは含めず、追加位置復元code/scale削除はmerge前残件。
