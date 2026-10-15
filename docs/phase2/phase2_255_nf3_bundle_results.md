# Issue #255：nf3四方式の束回転比較

## 実行契約と現在地

固定commit `0bed3d348fc74159ca218b2e292ace5a1db8fd19`。nf03 slots024、A/Bは追加frame位置復元scale1.25、C/Dは1.0（追加分なし）。A/Cはmotor横方向補正OFF、B/DはON。全方式とも解析torsion・energy-gradient frame反作用・局所one-ring実全vector反作用を共通使用する。排除体積・初期配置・hook剛性・長さ・T・dt・seedは固定。

ユーザーはcs10操作と実行を包括許可済み。queue #24の1τをレビュー通過し、queue #25の1sを2026-10-15 16:25:58 JSTに開始した。実効4 workers。1τ実測からの1s wall timeは約4時間02分、完了予定は20:30 JST頃＋負荷・変形依存。**現在は実行中で、goalは未達。2sとcanonical採択は行っていない。**

## 1τ short screen

#24は16:14:18～16:24:25 JST、全4条件各10,000 step・0.04s完走、集約成功。portable49ファイルをローカルへ同期、リモートとファイル集合・SHA-256を全件照合した。run_summaryを先に読み、共通development_evaluationで4/4 PASSを確認。全step support=5、solver success=1、fallback=0。固定camera3D/2D replay各1本（4方式panel、各41フレーム）と終端形状・時系列をレビューした。

| 方式 | QC | hook長相対誤差 最大 | motor torque残差比 最大 | 全系総トルク 最大 [N m] | 0.04sまでの菌体相対積算回転 [turn] |
|---|---|---:|---:|---:|---:|
| A | PASS | .03498 | 5.54e-16 | 4.69e-34 | .07634 |
| B | PASS | .03148 | 3.88e-16 | 2.51e-34 | .07666 |
| C | PASS | .06946 | 4.96e-16 | 4.22e-34 | .11464 |
| D | PASS | .06606 | 2.65e-16 | 2.10e-34 | .11657 |

回転計測valid fractionは全方式1.0。総力残差も全方式最大5.50e-28 N以下。motor横方向補正OFFでも、実トルクへの菌体側反作用により全系の釣り合いは成立している。

初期配置は120°/240°回転に対してRMS誤差4e-22m未満のC3対称。3本の積算回転が同一であることだけを自発的同期の証拠にしない。追加frame位置復元を外すと初期回転は増加するが、0.04sでは1回転未満で、束化・定常回転の採択判断には足りない。QC通過と急激な破綻がないことを根拠に全4方式を1sへ進めた。

## ローカル成果物

root：`outputs/2026-10-15/161410/issue255_nf3_bundle_execution/`

- `short_source/`：synced campaign、run_summary、state archive、rotation CSV
- `short_sync_verification.json`：49ファイルのSHA-256一致
- `short_evaluation/`：共通QC、初期配置、固定camera3D/2D replay
- `short_bundle_rotation/`：束の広がり・軸角・各べん毛回転の時系列/窓統計/図
- `short_stage_review.json`：short screenのlocal review PASSと全4方式の1s進行判断
- `manifest.json`：固定commit、queue、時刻、remote output root、実行許可

1s完了後に同じrootへ成果物を同期・照合し、全時間QC、束形成・軸整列・持続回転・駆動力・成分別釣り合いを比較して追記する。未完了段階で束回転成立を主張しない。
