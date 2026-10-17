# Issue #263：初期束配置と短期screen

## 比較契約

`conf/phase2_multi_run/2010_hex_nf3_bundle_initial_screen_issue263.yaml`の3条件を比較する。`nf03__slots024__C`は#255四方式比較のCと同じ初期座標・物理設定。`bundle`はopt-inの後方コンパクト配置、`bundle_phase5`は同配置のflag 1だけを**根元固定で菌体後方軸まわりに5°剛体回転**させた初期方位摂動である。設定名に`phase`が含まれるが、motor位相やらせん自身の軸まわりの位相ではない。3条件はnf3、slots024、T=2.5e-20 N m/本、dt_star=1e-4、Brownian/switching OFF、菌体–べん毛排除OFF、べん毛同士反発ONを共有する。

初期配置の手順・制限はADR 0027。preflightでは束半径・後方軸角・べん毛間軸角、hook長・角度・力、最小ビーズ間距離を記録する。非重複・hook長が成立しなければsimulationを始めない。t=0の幾何成立を遊泳の成立と扱わない。

## 実行契約と境界

Issue targetは`cs10_user_run`、独立3条件、parallel jobは`conf/phase2_parallel/issue263_initial_bundle/screen_job.yaml`。`cs10_qualified`、`max_workers: auto`、条件別出力を使用する。Macの長時間simulationは30分超と見積もり、実行しない。短期screenの実測後に長時間条件の所要時間を更新する。

Macではtargeted test、初期配置preview、`run_parallel.py --config conf/phase2_parallel/issue263_initial_bundle/screen_job.yaml --dry-run`を確認した。ユーザーは2026-10-17に本件のcs10接続・実行予約等の操作を明示許可した。過去の#255予約・許可は使用していない。

ローカルdry-run planは`outputs/2026-10-17/193744/dry_run_plan.json`、cs10でのdry-run planは同階層の`cs10_dry_run_plan.json`に保存した。実効3 workers、数値ライブラリthread各1、3条件の出力分離を確認した。cs10予約#26は固定commit `526a7f03ea276a9637715d060b30546a9f16df44`で開始し、既存dispatcherで実行した。別のdispatcherやtmuxは起動していない。

screen成果物は共通`development_evaluation`のfinite/body/hook length/flag/motor QCを先に確認し、固定cameraの3D／2D replayと`bundle_rotation`の時系列・周期診断を確認する。Cとの差分が初期幾何だけであることを保持する。QC未達なら幾何起因を切り分け、QC通過でも遊泳未達なら#255のA/B/C/D保存済み証拠と照合してCの物理前提を一項目ずつ再検討する。長時間stageの条件・必要回転周期数・挙動判定はscreenレビュー後に確定する。

## ローカル初期幾何確認

初期の探索では表面内向きhookを持つ配置が現れたため、候補をADR 0027の外向き・接線方向に限定した。制限後の画像とpreflightは`outputs/2026-10-17/193744/issue263_initial_geometry_preview/manifest.json`に保存した。初期束半径はC 1.041 µm、bundle 0.813 µm、phase5 0.807 µm。最小非付着ビーズ間距離は順に0.440、0.250、0.250 µmで、直径0.200 µmを超えた。hook方向は菌体外向き径方向から90°、hook長誤差は最大5.3e-23 m、hook初期力normは最大4.0e-28 N。これは静的確認であり、動的な成立は下記screenで別途評価した。

## 1τ短期screen：cs10予約#26

2026-10-17に3条件とも10,000 step・0.04 sを完走し、aggregate 3/3 completed、queue exit code 0。各`run_summary.json`を先に確認し、finite/body/shapeのオンラインゲートに失敗はなかった。ローカル`development_evaluation`はfinite/body/hook length/flag/motorを含む共通QCで**3/3 PASS**。固定cameraの3D／2D replay各1本（3条件panel、41フレーム）を生成し、終端3D像を確認した。共通evaluatorは固定値の`n_flagella`とmotor torqueをrun manifestから読めるよう補修した。QC閾値や物理設定は変えていない。

| 条件 | 束半径 中央値 [µm] | べん毛間軸角 中央値 [deg] | 速度 中央値 [µm/s] | 菌体相対積算回転 [turn] |
|---|---:|---:|---:|---:|
| C | 1.0133 | 2.650 | 0.492 | 0.11464 |
| compact rear | 0.7919 | 1.720 | 0.628 | 0.13516 |
| compact rear + flag 1初期方位5° | 0.7856 | 1.734 | 0.637 | 0.13478 |

束半径は全条件で初期より縮小し、終端は順に0.9932、0.7787、0.7724 µm。軸角は終端4.667、2.903、2.928°。初期方位摂動条件のflag 1対flag 0の積算回転差変動は0.0050 rad、回転診断のvalid fractionは全条件1.0。共通QCを通った保存状態のみを`bundle_rotation`へ渡した。完全な回転周期は全条件0件で、速度CVは0.354、0.413、0.405。これらは初期過渡を含む短い窓の診断であり、密な束の維持、位相同期、定常回転、持続推進を判定できない。

#255の保存済みA/B/C/D 1s比較では、frame追加位置復元を外したC/Dの束幅はA/Bより狭かった一方、Cでも2周期を通じて束幅が縮小し、収束は確認されなかった。今回のCの1τ回転0.11464 turnは#255短期Cと一致する。新配置の1τで束幅と速度が改善しても、Cのframe・motor前提の採択や長時間挙動の改善を意味しない。次はユーザーレビュー後、必要な完全回転周期数と長時間条件を決める。物理変更が必要ならframe復元・motor補正などを一項目ずつ対照化する。

成果物は`outputs/2026-10-17/195417/issue263_short_screen/`に保存。`source/`のportable 38ファイルをcs10とファイル集合・SHA-256全件照合した。`manifest.json`にqueue、固定commit、主要SHA-256、QC、replayを記録し、`evaluation/`と`bundle_rotation/`を同じrootに保持した。operational logは同期していない。

## 次タスク：nf=3～6を各1配置で比較（実行前契約案）

ユーザーの指摘により、次の比較から初期方位5°摂動を外し、nf=3～6の各本数で1配置ずつを試す。n3の旧3条件screenは履歴証拠として残すが、次stageの3 armをそのまま2秒へ延長しない。各本数の配置は、nf3 slots024、nf4 slots0134、nf5 slots01234、nf6 slots012345とし、全flagの初期方位offsetを0°に固定する。いずれもCのframe・motor・torsion等の設定、各flagのトルク、時間刻み、seedを共有する。本数と付着slotの組合せは同時に変わるため、1条件ずつの結果から本数単独の因果効果や成立範囲を断定しない。nf4～6は探索的診断であり、training scopeやcanonical採択の変更ではない。

| 本数 | 付着slots | 初期束半径 [µm] | 最小非付着ビーズ間距離 [µm] |
|---|---|---:|---:|
| 3 | 024 | 0.813 | 0.250 |
| 4 | 0134 | 0.813 | 0.250 |
| 5 | 01234 | 0.792 | 0.250 |
| 6 | 012345 | 0.813 | 0.250 |

上表はsimulationなしのgeometry preflight。ビーズ直径0.200 µmを超え、hook長誤差は1e-15 m未満、初期hook力normは1e-18 N未満だった。**静的成立は動的QCや束維持を保証しない。** n4～6は未screenの新条件なので、`conf/phase2_multi_run/2010_hex_nf3to6_bundle_1tau_issue263.yaml`の4条件（各0.04秒・10,000 step）を共通short-screen QCと固定camera replayに通す。parallel jobは`conf/phase2_parallel/issue263_initial_bundle/count_screen_job.yaml`。n3も同じjobに含めて同一commit・環境で比較し、旧#26のn3束条件との再現性を確認する。screenのQC未達条件は長時間へ進めず、Cおよび同一jobの他本数との対比較で原因を切り分ける。Macのsimulationなしdry-runは`outputs/2026-10-17/223411/issue263_count_plan/`に記録し、4条件の出力分離・4 workers・数値thread各1を確認した。

4条件のQCとreplayをレビューした後、通過条件だけを初期状態から2秒（50τ、500,000 step/条件）観測する候補を`conf/phase2_multi_run/2010_hex_nf3to6_bundle_2s_issue263.yaml`と`conf/phase2_parallel/issue263_initial_bundle/long_2s_job.yaml`に用意した。1τからのresumeではない。現時点のjobは4条件を列挙しているが、screen未達条件が出たら実行前に範囲・configを改め、予約を自動継続しない。`cs10_qualified`、autoの実効4 workers、各数値ライブラリthread 1、条件別出力、全条件geometry preflightを用いる。Macでは長時間simulationしない。旧nf3の実測外挿では2秒に8～11時間程度だが、nf4～6の計算負荷は未測定なので、screen実測を使って予約前にwall timeを更新する。

長時間archiveは各`run_summary.json`を先に読み、共通`development_evaluation`のfinite/body/hook length/flag/motor QCを判定し、PASSの保存状態のみを`bundle_rotation`へ渡す。各本数を固定camera 3D／2D replayで確認する。flag 0の菌体相対best-fit axial spinの完全周期を基準に、束半径、べん毛間軸角、輪郭長75%位置と先端の広がり、周期平均速度を同位相で比較する。これは厳密なmaterial twistや位相的な絡み合いの直接測定ではない。n6では対称性だけによる回転一致を同期の証拠としない。

**観測量の十分性**は本数ごとに完全周期5件以上を目標とし、先頭2周期を緩和として分離、末尾3周期の形状変化を記録する。#255のCは最初の2周期でも束半径が約9.7%縮小した。末尾の束半径または75%／先端広がりが周期間で5%を超えて変化した場合は、延長レビューの旗を立てる。これは新しい物理QCのFAIL閾値でも自動採択基準でもない。5周期未満・QC未達・継続的な形状変化があれば、成立／不成立を断定せず、観測延長または物理仮説の見直しをユーザーと判断する。

旧nf3・3条件の2秒案（commit `027403b`でdry-run済み）は**未予約・未実行で置き換え**る。旧planは`outputs/2026-10-17/222223/issue263_long_2s_plan/`に履歴として保存する。新しい4条件の固定commit・Mac/cs10 dry-run、実行範囲・所要時間を確認する前に予約・開始しない。新規hook・motor・相互作用の変更やcanonical freezeは含めない。
