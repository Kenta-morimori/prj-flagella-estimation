# Issue #263：nf3の5仮説・1.5秒比較

## 許可と実行範囲

2026-10-27、ユーザーは5条件の比較設計を概ね了承し、cs10操作・予約を許可した。直後の確認で実行長は「1.5h」ではなく **1.5秒**と確定した。物理モデルの採択、merge、追加の条件探索は含めない。

Issue targetは`cs10_user_run` / `execution:cs10`。PR #264のfeature branchで実装し、primary checkoutの別Issue作業は変更しない。Macでは短いtargeted test、full regression、simulationなしpreflight/preview/dry-runのみ。新しい候補は先行1τのstrict QCと固定camera replayをレビューしてから本比較へ進む。FAILなら報告して停止し、係数・閾値・刻み変更や自動再実行をしない。

## 5条件と共通値

| ID suffix | 単独の変更 |
| --- | --- |
| `physical_baseline` | B。修正済み線分距離・従来compact・分布motor・指数cutoff反発 |
| `local_motor` | M。`motor.segment_force_placement: root_three`。同じ形状・同じ拡散分配状態でBの実全トルクvectorを根元3 beadへ移す。軸・全ベクトル菌体反作用は一致。退化時は停止し支持点を変更しない |
| `continuous_repulsion` | C。`potentials.spring_spring_repulsion.force_curve: continuous_cutoff`。d→0+での力を維持し、cutoffで連続して0となる力曲線へ変更 |
| `root_relocated` | G。`flagella.initial_root_tangent_angle_deg: 39.276`。後方軸からaxis×outwardの接平面方向へのhook方向角。各らせんを剛体並進し、形状・位相・hook長を保持。0は追加移動なし |
| `body_repulsion_on` | R。Bから菌体–べん毛反発だけON |

全IDのprefixは`nf03__slots024__`。nf3/slots024、T/reference T=2.5e−20 N m/本、内部dt=4 µs、粘性1 mPa s、菌体30 bead・各鎖11 bead、hook長0.25 µm、物性・解析torsion・Brownian OFF・motor切替OFF・保存1 msを固定する。べん毛間反発は常時ON。B/M/C/Gの菌体反発はOFF。追加位置復元scaleは全条件1.0、腕重みと横トルク補正は使わない。

M/C/Gはopt-inで、既定profileの値を変更しない。Mは拡散状態をBと同じ方法で更新し、分布motorが作る全トルクを局所最小ノルムの力対へ置き換える。菌体反作用は元の全トルクから計算するため、同一状態でBと完全一致する。これは力の作用位置の機序比較であり、根元接線軸の別モデルへの変更ではない。

Cの式・Gの静的設計と検証限界は[比較計画](phase2_263_nf3_physical_comparison_plan.md)に記録する。Cはcutoff内の力曲線も変えるので不連続だけの因果対照ではない。Gのbeadクリアランス保持は連続面/segmentの非貫通保証ではない。いずれも改善を前提としない。

## checked-in実行設定

| 段階 | campaign | parallel job | 条件別step |
| --- | --- | --- | ---: |
| 1τ=0.04秒 | `conf/phase2_multi_run/2010_hex_nf3_physics_1tau_issue263.yaml` | `conf/phase2_parallel/issue263_initial_bundle/nf3_physics_1tau_job.yaml` | 10000 |
| 37.5τ=1.5秒 | `conf/phase2_multi_run/2010_hex_nf3_physics_1p5s_issue263.yaml` | `conf/phase2_parallel/issue263_initial_bundle/nf3_physics_1p5s_job.yaml` | 375000 |

本比較計1875000 step。runtimeはtau=0.039999999999999994秒、duration=1.4999999999999998秒、各375000 stepとなる。過去2秒の500001 stepとは異なり、今回の設定で余分なstepは予定されない。

両jobとも`max_workers: auto` / `cs10_qualified`、実効5 workers、数値thread各1、条件別の子run output namespaceを使う。先行screenをレビューするまで本比較を開始しない。queueにこのゲートを自動判定する依存機能はないため、本比較を先に無条件enqueueしない。

過去nf3 ONの4.718 step/sから本比較の最遅参考wallは約22.1時間（queue待機・準備・先行screen・解析を除く）。距離修正後や新候補の速度は未測定なので、先行screenの実測で更新する。

## 検証と証拠

- `tests/test_nf3_physical_comparison.py`：5条件が単独差分であること、1τ→1.5秒の時間だけの延長、Mの力/トルク保存と菌体反作用一致・退化時停止、Cのポテンシャル勾配/連続性、Gの形状/位相/hook/対称性、各engineの短いintegrationを検証する。
- full regression：1094件PASS（209.58秒）。ruff format/checkとdiff checkもPASS。
- targeted suite：物理比較・potential・motor・既存Issue263初期配置・固定bundleの84件PASS。
- 距離修正実行commit `0017fbbd922997f32dbad4b40aa5af03da974c92`とのBの8内部step比較：座標と6力成分、計56配列がbitwise一致。標準の分布motorと指数反発を保持している範囲の回帰確認である。
- 共通geometry preflightはcompactでも初期弾性/反発力を記録する。Gの軸ずれ0.260039→0.156806 µm、各鎖並進0.168038 µm、hook長0.25 µm・最小bead距離0.25 µmを確認する。Rの初期反発は隠さず記録する。
- ローカル準備証拠：`outputs/2026-10-27/174300/issue263_nf3_physics_preparation/`。両dry-run JSON、初期preview、baseline regression、full pytest log、manifestを保全する。固定commitとreservationはenqueue後のmanifest/PRへ記録し、記録のためだけに実行commitを置換しない。

## 完走後

固定commit・5ID・job/campaign SHAと予定stepを照合し、portable outputsだけをJST日時のローカルoutputsへ同期する。ファイル集合/件数/全SHAを照合し、各`run_summary.json`を先に読む。共通`development_evaluation` strict QC・heatmap・固定camera同一視野3D/2D replayを生成する。

QC-PASS archiveで各鎖材料点spin、菌体材料点自転、束幅/軸整列、符号付き前進/停止/後退、根元往復と姿勢、接触、前進/回転、motor仕事を同時刻・可能な同位相で比較する。初期0〜0.2秒、既往減速窓0.4〜0.8秒、後半0.8〜1.5秒と初期2回転を除く完全周期を分ける。周期数不足やQC FAILは未達として残し、時間だけで定常性を採択しない。

Issueには重要結果・未達・次の判断を3〜5項目で既存要約コメントへ反映し、詳細は本runbook/結果へリンクする。PR #264・local review_resultも更新する。新たな物理変更や次段階の長時間jobは別判断である。
