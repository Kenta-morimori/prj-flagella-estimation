# Issue #263：初期束配置と短期screen

## 比較契約

`conf/phase2_multi_run/2010_hex_nf3_bundle_initial_screen_issue263.yaml`の3条件を比較する。`nf03__slots024__C`は#255四方式比較のCと同じ初期座標・物理設定。`bundle`はopt-inの後方コンパクト配置、`bundle_phase5`は同配置のflag 1だけを後方軸まわりに5°回す。3条件はnf3、slots024、T=2.5e-20 N m/本、dt_star=1e-4、Brownian/switching OFF、菌体–べん毛排除OFF、べん毛同士反発ONを共有する。

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
| compact rear + flag 1位相5° | 0.7856 | 1.734 | 0.637 | 0.13478 |

束半径は全条件で初期より縮小し、終端は順に0.9932、0.7787、0.7724 µm。軸角は終端4.667、2.903、2.928°。位相摂動条件のflag 1対flag 0の積算位相差変動は0.0050 rad、回転診断のvalid fractionは全条件1.0。共通QCを通った保存状態のみを`bundle_rotation`へ渡した。完全な回転周期は全条件0件で、速度CVは0.354、0.413、0.405。これらは初期過渡を含む短い窓の診断であり、密な束の維持、位相同期、定常回転、持続推進を判定できない。

#255の保存済みA/B/C/D 1s比較では、frame追加位置復元を外したC/Dの束幅はA/Bより狭かった一方、Cでも2周期を通じて束幅が縮小し、収束は確認されなかった。今回のCの1τ回転0.11464 turnは#255短期Cと一致する。新配置の1τで束幅と速度が改善しても、Cのframe・motor前提の採択や長時間挙動の改善を意味しない。次はユーザーレビュー後、必要な完全回転周期数と長時間条件を決める。物理変更が必要ならframe復元・motor補正などを一項目ずつ対照化する。

成果物は`outputs/2026-10-17/195417/issue263_short_screen/`に保存。`source/`のportable 38ファイルをcs10とファイル集合・SHA-256全件照合した。`manifest.json`にqueue、固定commit、主要SHA-256、QC、replayを記録し、`evaluation/`と`bundle_rotation/`を同じrootに保持した。operational logは同期していない。

## 次タスク：2秒の長時間観測チェックポイント（実行前契約案）

1τで完全周期は0件。#255のCは1秒で2.519周、最初の2周期で同位相の束半径が約9.7%縮小した。巻き付き・束形状の緩和と、その後の維持を分けるため、`conf/phase2_multi_run/2010_hex_nf3_bundle_2s_issue263.yaml`でC／compact rear／flag 1位相5°の**同じ3条件を初期状態から2秒**観測する。1τからのresumeではない。3条件間の物理設定・seed・時間刻み・出力間隔は短期screenと同一で、時間指定だけ2秒（50τ、500,000 step/条件）にする。#255の旧C 1秒archiveは履歴比較に使うが、2秒の直接対照を代替しない。

`conf/phase2_parallel/issue263_initial_bundle/long_2s_job.yaml`は`cs10_qualified`、autoの実効3 workers、数値ライブラリ各1 thread、条件別出力と全条件geometry preflightを使う。Macではsimulationしない。#26の1τ約10分および#25の1秒約4時間45分からの単純外挿では、cs10 wall timeは概ね8～11時間を見込む。共有機負荷で変動するため保証値ではない。2秒は**中間判定点**であり、ここで定常成立を自動判定しない。

Mac上のsimulationなしdry-runは`outputs/2026-10-17/222223/issue263_long_2s_plan/dry_run_plan.json`に保存した。3条件の出力分離、preflight、実効3 workers、数値thread各1を確認済み。固定commitでのcs10側dry-runは実行開始前に照合する。

判定の順序と記録を実行前に固定する。

1. 各`run_summary.json`を先に読み、3/3完走と共通`development_evaluation`のfinite/body/hook length/flag/motor QCを判定する。QC未達条件のarchiveを束診断へ流さず、Cとの対比較で初期幾何起因を調べる。
2. QC通過archiveから固定cameraの3D／2D replay、`bundle_rotation`時系列・完全周期CSVを生成する。完全周期はflag 0の菌体相対best-fit axial spinで定義し、位相欠測・逆転・保存不足を「判定不能」として残す。これは厳密なmaterial twistでも、3本の自発同期の直接証拠でもない。
3. 各周期の同位相境界で束半径、べん毛間軸角、輪郭長75%位置・先端の横方向広がりを比較する。周期平均の回転周期・速度、flag間の積算位相差、推進方向も併記する。輪郭広がりは接触や位相的な絡み合いを直接証明しないので、replayを併せて確認する。
4. **観測量の十分性**として各条件で完全周期5件以上を目標にし、先頭2周期を初期緩和として分離して末尾3周期を比較する。#255のCでは2周期でも形状が変化し続けたため、この数を採る。5周期未満、QC失敗、または末尾3周期でも明瞭な一方向の形状変化が続けば、成立／不成立を断定せず観測延長または仮説修正をレビューする。

最後の3周期では、同位相の束半径と75%／先端広がりの周期間変化率を数値で示す。5%を超える変化は#255の約9.7%縮小を見逃さないための**延長レビュー旗**とする。新しい物理QCのFAIL閾値やモデル採択基準にはしない。束軸・べん毛間軸角、周期、速度、位相差は各周期の値と傾向を提示し、許容変動の採択はデータとreplayのレビューで判断する。2秒で十分な周期が取れない場合、3～4秒への延長候補を別の固定commit契約として準備し、予約を自動追加しない。

この節は実行条件をレビュー可能にする準備であり、長時間jobの予約・開始を記録するものではない。長時間開始前に、固定commit、job範囲、dry-run plan、所要時間を示してユーザー判断を受ける。新規hook・motor・相互作用の変更、canonical採択、本数展開は含めない。
