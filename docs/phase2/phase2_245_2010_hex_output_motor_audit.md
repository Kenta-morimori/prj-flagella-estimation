# Issue #245: 2010 hex 2秒archiveの出力・motor torque診断（2026-09-28）

対象は`2010_hex_project`の同期済み13条件・2秒completed archive。simulation、cs10操作、特徴量解析、採択判定は実施していない。現行`development_evaluation.long_duration`のstrict判定は**13/13 FAILのまま**であり、以下はdiagnostic-onlyの結果である。

## 出力

- 3D replayの共通既定値をbody-follow OFF（固定）にした。明示的な`--camera-3d follow`とprofile側の明示ONは維持する。13条件をn別6群、各41 frameで保存済みarchiveから再生成し、全manifestで`camera_3d=fixed`、`view_range_mode=campaign-envelope`を確認した。旧#245 replayも同じ明示設定だったため、6群すべての最終PNGは旧版とSHA-256が完全一致し、構図の変化はない。2D replayは再生成していない。3Dの番号・色・QC表示も変更していない。
- 新しい共通の付着slot平面図は13条件を一枚に配置する。初期菌体の後方（−x側）から+x向きに中央六角環を見る。横軸+y、縦軸+z、0始まりのslot番号、3D replayと同じF番号・色を表示し、実archiveの初期body bead位置と保存済み`geometry.actual.attachment_topology`を照合した。n=6は`full_ring_rotation_equivalent`である。

ローカル成果物：`outputs/2026-09-28/215511/issue245_output_review/`（`evaluation/attachment_slots/attachment_slots_all_conditions.png`、同`manifest.json`、`replay3d/nf01`〜`nf06`のMP4・最終PNG・manifest・run.log）。入力archiveのSHA-256は統合summaryおよびslot図manifestに記録した。

### 固定世界座標での菌体移動を見せる追加replay（2026-09-29）

`conf/sim_swim_2010_hex.yaml`には`render.follow_camera_3d`の指定がない。旧#245 replayも`--camera-3d fixed`と`campaign-envelope`を明示しており、全13 archiveのビーズ位置から一度決めた共通カメラ中心・画角を各frameに固定していた。したがって旧動画が菌体を追従していたわけではない。

ただし旧画角幅は約7.33 µmで、菌体の初期→最終中心移動0.18〜1.36 µm（中央値0.61 µm）が視覚的に小さく見える。そこで**世界原点(0,0,0) µm固定、各軸±2.8 µm固定**の追加3D replayを、同じ13 completed archiveから生成した。初期から最終まで全body beadが画角内にあることを検証した。遠位のflagellaは画角外に出る可能性があるため、全形状観察には元のcampaign-envelope版を使う。追加動画もdiagnostic-onlyであり、strict判定・3D色・番号・QC表示・2D replayは変更しない。

追加成果物：`outputs/2026-09-29/095317/issue245_fixed_world_body_motion/replay3d/nf01`〜`nf06`のMP4・最終PNG・manifest・run.log。`view_range_mode=explicit-fixed`は3D fixed camera専用のgeneric replay CLI optionで、追従カメラや既存profileの明示設定は変更しない。

## Motor torqueの同一状態診断

現行の`root_torque_segment_couples`とnominal local-twist重みを再構成し、各条件の初期状態、最初の**観測済み**0.02超過、診断sample中の最大残差、最終観測状態を選んだ。初期状態には記録済みmotor診断sampleがないため、反実仮想計算は保留した。観測sampleとの照合はbody/flagの3D torque vectorの相対誤差≤0.005、torque/force残差比の絶対誤差≤0.005を条件とした。

観測状態の選択39件中16件で再計算値が記録値に一致した（最初の超過13/13、最大残差1/13、最終2/13）。一致した状態だけで全bodyに全ベクトル反作用を与える同一geometry・同一重みの比較を実施した。現行の再計算torque残差比は0.0208〜0.2925、全ベクトル反作用では最大`2.54e-16`だった。**これは力の釣り合いを代数的に改善できるという診断であり、変更後の運動、束軸安定性、2秒の実行成功を示さない。** 残り23件は記録済みvectorとの不一致により反実仮想を出していない。保存archiveの観測時刻とstep内の`positions_before_m`には差があり得るが、不一致の原因は本診断のみでは確定できない。

ローカル成果物：`outputs/2026-09-28/215511/issue245_output_review/motor_torque_audit/manifest.json`と`run.log`。manifestは各Fの軸方向・横方向torque、菌体反作用、全体の力・torque残差、記録値照合、比較の有無、入力archiveとdiagnostic samplesのSHA-256を保持する。

### 同一stepの記録値による時系列再解析（2026-09-29）

archiveと診断値の時点差を避けるため、13条件の`diagnostic_samples.csv`に**同じstepで記録された**菌体側・べん毛側の3D motor torque vectorを直接合算した。各condition 201 sample、計2,613 sampleであり、10 ms間隔の観測値である。`run_summary.json`の全step最大値と照合し、completed archive・`run_summary.json`・`performance.json`のSHA-256を統合評価表に照合した。入力診断CSVと出力図表のSHA-256は追加manifestに記録した。

全13条件で記録sampleが`0.02`を超え、最初の**sampled**超過は0.01–0.09 sである。2,613 sample中2,490 sampleが閾値超過し、sampled最大残差比は条件別に0.1705–0.3049、全step最大値は0.1990–0.3661であった。sampled最大は全step最大の81.5–96.5%に留まるため、この図から真の初回超過stepや最大stepを推定しない。記録されたmotor合力残差比のsample最大は`1.33e-16`で、力の不釣り合いではない。

現行のtorque residual ratioは`||body torque + flag torque|| / Σ_i ||(r_i-r_mean)×F_motor,i||`であり、「目標motor torqueの何%か」ではない。追加CSVでは別途`||body+flag|| / (n×|T|)`も出し、sampled最大時点で条件別に1.54–7.19となった。正規化の違いを保ったまま両者を比較し、閾値`0.02`の物理的意味づけは別途判断する。

最初のsampled超過時点は13/13条件で既存archive再構成と記録値が照合可能だった。計42本のべん毛について、各軸方向の`flag + body`残差はnominal torque比で最大`1.25e-15`、軸直交成分は0.275–1.561倍だった。これは**現行方式が軸方向のみを菌体に反作用として与え、実際に生成された横方向トルクを相殺していない**ことを支持する。保存状態での全ベクトル反作用が代数的に残差を消す結果とも整合する。ただし他の選択観測23件ではarchive再構成が記録値と一致していない。実装上、diagnostic torqueはstepの`positions_before_m`で計算され、compact archiveはstep後の観測状態を保存し、archive境界では補間する。これが不一致の有力要因だが、23件それぞれの原因と大きさは同一stepのforce snapshotがないため確定できない。

追加成果物：`outputs/2026-09-29/101619/issue245_motor_torque_recorded_analysis/`の`recorded_torque_timeseries.png`、時系列・condition summary CSV、照合済み初回超過のF別torque CSV、manifest、run.log。これはdiagnostic-onlyであり、現行strict FAIL 13/13と閾値`0.02`を変更しない。`body_reaction_full_vector=true`は全body beadへ反作用を分配するproject実装であり、局所hook反作用の物理的適切さ、修正後の軌道・束軸・長時間安定性は別途短時間比較で検証する必要がある。

[Watari & Larson (2010)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2800969/)はhook近傍のtorque釣り合いとredirecting torqueを記述する。一方、このprojectの`body_reaction_full_vector=true`は反作用を**全body bead**に分配する実装上の反実仮想であり、論文の局所hookモデルと同一ではない。現行strict閾値0.02と力モデルは変更していない。motor residualの物理的妥当性、全ベクトル反作用での再実行、束軸時間変化の定量化は後続タスクとする。
