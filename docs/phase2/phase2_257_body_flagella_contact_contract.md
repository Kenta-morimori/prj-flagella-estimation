# Issue #257: body--flagella 排除相互作用比較契約

## 目的と範囲

`2010_hex_project`を主対象、`2010 project`を補足対象として、body-only segmentとflagellumを含むsegmentの排除相互作用を無効化した場合に、多べん毛条件の一時停止が減るかを診断する。flagella--flagella反発は維持し、body--bodyは既存の拘束ベース挙動を維持する。この比較はcanonical model、dataset、既存strict QC閾値、#255のmotor反作用方針を変更しない。

## 実行前のtriage

Issue #257は`execution:triage`である。short screenを含む実行、cs10への接続、reservation操作は、実行先・実測wall time・並列job・dry-runを確定し、ユーザーが当該操作を明示承認するまで行わない。

## 比較設定

- hex short screen: `2010_hex_project_body_flagella_contact_screen_issue257.yaml`。#245の13 attachment配置をON/OFFで比較する26条件。
- hex 2 s: `2010_hex_project_body_flagella_contact_2s_issue257.yaml`。short screenをreviewし、ON/OFF両armがfinite/body/hook length/flag/motor gateを満たすpairだけを対象にする。
- project補足: `2010_project_body_flagella_contact_{screen,2s}_issue257.yaml`。`n=4,5,6`をON/OFFで比較する各6条件。
- 全armでtorque、Δt、seed、archive、Brownian/switchingを固定する。motor residualは#255との並行診断として記録し、採択根拠にはしない。

## 2秒診断

40 ms窓でbody速度とbody-rollが各condition中央値の10%未満となる同時区間をstall候補とする。body長軸角度は姿勢揺らぎ・方向転換の補助時系列とし、固定camera 3D replayで束化中の連続回転と整合するかを確認する。排除OFFでstall頻度または総時間が減少しても、仮説を支持する診断結果に留め、model採択・strict QC PASSとは解釈しない。
