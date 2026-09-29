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

[Watari & Larson (2010)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2800969/)はhook近傍のtorque釣り合いとredirecting torqueを記述する。一方、このprojectの`body_reaction_full_vector=true`は反作用を**全body bead**に分配する実装上の反実仮想であり、論文の局所hookモデルと同一ではない。現行strict閾値0.02と力モデルは変更していない。motor residualの物理的妥当性、全ベクトル反作用での再実行、束軸時間変化の定量化は後続タスクとする。
