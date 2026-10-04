# Issue #257: body--flagella 排除相互作用比較契約

## 目的と範囲

`2010_hex_project`を主対象、`2010 project`を補足対象として、body-only segmentとflagellumを含むsegmentの排除相互作用を無効化した場合に、多べん毛条件の一時停止が減るかを診断する。flagella--flagella反発は維持し、body--bodyは既存の拘束ベース挙動を維持する。この比較はcanonical model、dataset、既存strict QC閾値、#255のmotor反作用方針を変更しない。

## 実行前のtriage

Issue #257は`execution:cs10`のUser-run campaignである。short screenを含む実行、cs10への接続、reservation操作は、並列job・dry-runを確定し、ユーザーが当該操作を明示承認するまで行わない。

## 比較設定

- hex short screen: `2010_hex_project_body_flagella_contact_screen_issue257.yaml`。#245の13 attachment配置をON/OFFで比較する26条件。
- hex 2 s: `2010_hex_project_body_flagella_contact_2s_issue257.yaml`。short screenをreviewした後、#245のcompleted ON 13条件を再利用し、OFF 13条件だけを新規実行した。motor strict FAILが残るため、ユーザー承認により両armを診断専用で比較した。
- project補足: `2010_project_body_flagella_contact_{screen,2s}_issue257.yaml`。`n=4,5,6`をON/OFFで比較する各6条件。
- 全armでtorque、Δt、seed、archive、Brownian/switchingを固定する。motor residualは#255との並行診断として記録し、採択根拠にはしない。

## 2秒診断

40 ms窓でbody速度とbody-rollが各condition中央値の10%未満となる同時区間をstall候補とする。body長軸角度は姿勢揺らぎ・方向転換の補助時系列とし、固定camera 3D replayで束化中の連続回転と整合するかを確認する。motor residualはstrict FAILとして保持するが、ユーザーreview後の2秒診断を停止させない。排除OFFでstall頻度または総時間が減少しても、仮説を支持する診断結果に留め、model採択・strict QC PASSとは解釈しない。

## 2秒ON/OFF比較結果（2026-10-04）

#245のcompleted ON 13条件と#257のcompleted OFF 13条件を結合した共有workspace bundleは`outputs/2026-10-04/131754/issue257_hex_2s_on_off_visualization/`にある。13ペアすべてで上記stall定義の候補はON/OFFとも0窓・0秒だった。ゼロ同士の比較から排除OFFによるstall改善の有無は判定できず、本結果を仮説の支持・反証やmodel採択に用いない。nf=5のreplayには回転・並進・姿勢の乱れが見られるが、その機構は特定していない。

26条件すべてで`first_failure_category=hook`が最初の内部step（4 µs）に記録された。元の生記録はhook角度誤差30°超を拾うが、角度は本契約ではdiagnostic-onlyである。実行設定から初期モデルを再構成すると、t=0の角度誤差は全26条件で31.106〜32.383°であり、ON/OFFで同一だった。したがって初回記録は初期配置と角度基準の不整合に由来し、排除ON/OFFの差や最初の時間積分による発生ではない。既存archiveの形状・閾値は変更しない。

hook長の最大相対誤差は全条件で0.0509〜0.1007とstrict閾値1.0未満である。一方、motor torque residualは0.1007〜0.3661でstrict閾値0.02を全条件で超え、26/26 armがstrict FAILである。比較はdiagnostic-onlyのままとし、#255の解決根拠、canonical model、datasetには渡さない。projectの2秒比較は未実行であり、#257は本PRのmergeでは完了しない。

## pending初期hook中立配置

PR #259で既定OFFの`flagella.initial_hook_force_neutral`を追加し、`2010_hex_project_hook_neutral_screen_issue257.yaml`の26条件だけでONにする。各べん毛の位相・螺旋軸・内部bond/曲げ/ねじれを変えずに全体を平行移動し、hook長を固定して初期hook角を90°にする。90°はhookポテンシャルの適用分岐に入るが平衡角であるため、t=0の曲げ力は数値許容差内でゼロとなる。preflightでは角度90°±1e-6°、hook長、hook力、外向き、非付着body beadとflagellum beadおよびflagellum同士の中心距離≥bead直径を確認する。bead clearanceはsegment間の完全非接触を意味せず、後者は固定した形状・hook長の下では保証しない。

この変更は初期geometryのみの候補であり、simulationは未実行である。旧#245/#257の26 archiveは修正前モデルの診断結果として分離し、新候補へのQC PASSや長時間結果として再解釈しない。従来の30° raw指標、動作中のstrict QC閾値、motor residualの扱いは変えない。新short screenの実行には別途ユーザーの操作別承認が必要である。
