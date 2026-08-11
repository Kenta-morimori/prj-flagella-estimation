# Phase 2 Decision Record

このファイルは，Phase 2で採用・不採用・置換・保留した判断を，後続開発で再利用できる形で記録する．

Issue単位の進捗台帳，branch一覧，acceptance criteria一覧，実行command一覧ではない．  
現在の作業状態は`phase2_current.md`，詳細な実装契約はconfig・schema・test・task-specific doc，設計理由はADRを正本とする．

## Status

- `adopted`: 現在採用している
- `rejected`: 評価したが採用しなかった
- `replaced`: 過去に採用したが，後続判断で置換した
- `diagnostic`: 診断結果として保持するが，canonical条件ではない
- `pending`: 評価またはユーザー判断が未完了

## Decision Index

| ID | Decision | Status |
|---|---|---|
| P2-D01 | 初期geometryの検証契約 | adopted |
| P2-D02 | body・hook・single flagellumの段階的stability gate | adopted |
| P2-D03 | 螺旋net回転の評価方法 | adopted |
| P2-D04 | root torqueのflagellum chainへの伝搬方式 | adopted |
| P2-D05 | 2010 projectのmotor mode・torque帯・local scale方針 | adopted |
| P2-D06 | 複数べん毛の軸整列・束化判定 | adopted |
| P2-D07 | attach-frame補強とbasal自由度 | adopted |
| P2-D08 | べん毛本数分析のfeature・dataset contract | adopted |
| P2-D09 | seed依存初期配置とcenter-priority policy | adopted |
| P2-D10 | multi-run・sweep・raw output・replay policy | adopted |
| P2-D11 | diagnostic dataset v0の位置づけ | diagnostic |
| P2-D12 | n≥4安定化とdataset v1の採択範囲 | adopted |
| P2-D13 | dataset version・revision・provenance | adopted |
| P2-D14 | Phase 3 handoffと2D識別性の初期判断 | adopted |
| P2-D15 | 2010 potential formulationとproject初期軸 | adopted |
| P2-D16 | model profileと時間schema | adopted |
| P2-D17 | 2015 refined geometryとpaper-inspired motor | adopted |
| P2-D18 | 2015 refined model Stage A採否 | pending |
| P2-D19 | dataset v2前のn=3長時間failure診断 | pending |
| P2-D20 | RUN–TUMBLEの段階実装 | pending |

---

## Geometry・stability・motor model

### P2-D01: 初期geometryを自動検証可能なcontractとして扱う

- **Status:** adopted
- **Background:** 参照論文に基づくbody・hook・flagellumの初期形状が，simulation開始前から妥当であることを定量確認する必要があった．初期geometryの誤差と時間発展後の破綻を分離できなければ，failure原因を正しく判断できない．
- **Change:** bond length，helix pitch・radius，bending，torsionなどのtarget・tolerance・pass/failを`initial_geometry_summary.json`へ出力し，初期stateとstep summaryの整合をtestで固定した．
- **Result:** 2010 geometryについて，論文由来値，project実装値，許容範囲を機械的に検証できるようになった．
- **Interpretation:** 初期geometry gateは，motor・stability・dataset評価より前に通過すべき基礎contractである．詳細なparameter値と許容値はconfig・builder・testを正本とし，本ファイルへ複製しない．
- **Decision:** 初期geometry検証をPhase 2の恒常的な入口gateとして維持する．
- **Evidence:** Task P2-2-001，`tests/test_simulation.py`，`tests/test_model_builder.py`．

### P2-D02: failureをbody・hook・flagellumの段階的gateで切り分ける

- **Status:** adopted
- **Background:** full swimmerを一度に評価すると，body変形，hook過伸長，flagellum形状破綻のどこがfirst failureか判別しにくかった．
- **Change:** body-only，body + minimal basal stub，body + hook + single full flagellumの順に複雑度を上げ，safe representativeとfirst-fail representativeを固定した．`shape_pass_nonbody`，body gate，`first_fail_category_nonbody`などを用いてfailure主体を区別した．
- **Comparison:** torqueなし／あり，低torque／高torque，minimal basal stub／single full flagellumを比較した．
- **Result:** body-onlyの安定帯，hook first-fail boundary，single flagellumでflag bond・bend・torsionが支配的になる条件を再現可能にした．
- **Interpretation:** full swimmerの破綻は，body・hook・flagellumを分けて診断すべきであり，最終行だけでなくfirst-fail categoryを用いる必要がある．
- **Decision:** 段階的stability gateとfirst-fail分類を維持する．過去の個別representativeはregression testまたはdiagnostic evidenceとして扱い，canonical model条件とは区別する．
- **Evidence:** Tasks P2-3-002，P2-4-003，P2-5-004，`tests/test_simulation.py`，`tests/test_motor_scale_sweep.py`，`tests/test_motor_forces.py`．

### P2-D03: 螺旋回転は瞬間rateではなくnet回転と方向一貫性で判定する

- **Status:** adopted
- **Background:** root azimuth由来の`flag_phase_rate_hz`や`median(abs(flag_helix_spin_rate_hz))`は，fit jitterや往復揺れを持続回転として誤判定する場合があった．
- **Change:** `flag_helix_spin_phase_deg`の累積差からnet revolutionを計算し，回転方向一貫性とshape gateを組み合わせた．
- **Comparison:** 旧triplet条件では，0.25 sでrootが約0.88回転してもhelix全体は約0.0014回転に留まり，helix/root比は約0.0016だった．瞬間rateが高くても持続的な螺旋回転ではなかった．
- **Result:** 形状維持と螺旋全体のnet回転を独立に評価できるようになった．単一べん毛代表gateとして，net 1回転以上，方向一貫性0.5以上を用いた．形状側ではhook，flag bond，bend，torsionの既存gateを併用した．
- **Interpretation:** root回転，body roll，螺旋fit位相，螺旋軸中心回転は同一ではない．回転判定には累積phaseと方向一貫性が必要である．
- **Decision:** net revolutionとdirection consistencyを螺旋回転判定の主指標とし，瞬間rateは補助診断とする．
- **Evidence:** Tasks P2-6-005，P2-6-006，`tests/test_helix_retention_gate.py`，`tests/test_run_state_fixed.py`．

### P2-D04: root torque伝搬には`root_torque_segment_couples`を用いる

- **Status:** adopted
- **Background:** bead-position-onlyのtriplet motorでは，root方位は変化しても，flagellum chain全体へ軸方向torqueが十分に伝わらなかった．material frame，segment twist，axial torque fluxに相当する明示的な状態量が不足していた．
- **Change:** `axial_torque_flux_probe`，`local_twist_transmission_probe`，全体軸投影方式を比較した後，segmentごとの局所force coupleでroot torqueを伝える`root_torque_segment_couples`を実装した．
- **Comparison:** `triplet`はroot付近の回転・変形に留まりやすく，`root_torque_axis_projection`は全体へtorqueが届く場合の比較・診断に有用だった．probe modesは仮説検証には有効だが正式modeにはしなかった．
- **Result:** 0.5 s代表条件でshape gateを保ちながらhelix net 1回転以上を確認した．既存torsion forceをOFFにするとshape gateが破綻した．
- **Interpretation:** 既存torsion potentialは螺旋形状維持，`root_torque_segment_couples`はroot torque伝搬という異なる役割を持つ．新方式は参照論文を完全再現するものではなく，project-specific extensionである．
- **Decision:** 2010 projectの正式modeとして`root_torque_segment_couples`を採用する．`triplet`と`root_torque_axis_projection`は比較用として残し，probe modesは正式実行modeから削除する．既存torsion forceは維持する．
- **Evidence:** Tasks P2-6-007，P2-6-008，ADRs 0002，0003，0004．

### P2-D05: 2010 projectでは標準torque帯とcanonical motor名を固定する

- **Status:** adopted
- **Background:** `root_torque_segment_couples`導入後に，単一べん毛の安定torque帯とlocal motor scaleの必要性を再評価する必要があった．また，過去のprobe名と正式mode名が混在していた．
- **Change:** torque sweepとlocal spring・bend・torsion・hook scaleのone-factor比較を行い，canonical mode名を整理した．
- **Result:** 単一べん毛，0.5 s，`dt_star=1.0e-4`では`2.0e-20`〜`3.0e-20 N m`が代表PASS帯で，`3.5e-20 N m`以上ではflag failureが支配的だった．各local motor scaleを2.0へ上げても高torque failureを救済しなかった．代表値として`2.5e-20 N m`を採用した．
- **Interpretation:** 高torque failureは単純な局所stiffness不足ではなく，motor local scaleの増強を標準解とすべきではない．
- **Decision:** motor local scaleは原則`1.0`をbaselineとし，非`1.0`は診断用project extensionとする．canonical名を`triplet`, `root_torque_segment_couples`, `root_torque_axis_projection`とし，旧名はdeprecated aliasとして扱う．
- **Evidence:** Tasks P2-6-009，P2-6-010，Issues #54，#88，ADR 0004．

### P2-D06: 複数べん毛の束化は，まず螺旋軸の安定整列として判定する

- **Status:** adopted
- **Background:** 近接接触だけを用いたbundle判定では，明示的に接触していなくても複数べん毛軸が揃う有効な遊泳候補を除外する可能性があった．
- **Change:** 各flagellumの第2ビーズ以降から螺旋中心軸を推定し，pair angle，平均軸からの偏差，rearward angle，alignment orderを記録した．
- **Comparison:** close pairやbundle participationだけの判定と，軸方向整列を比較した．
- **Result:** 後半80%で各flagellum軸が平均軸から15 deg以内に揃うことを主なaxis-aligned判定とした．第1ビーズはbasal geometryの影響が強いためaxis PCAから除外した．
- **Interpretation:** axis alignmentは近接束化の完全な代替ではないが，短時間simulationで後方整列候補を比較する再現可能な指標である．hook failureやshape failureは別gateとして併記する必要がある．
- **Decision:** axis alignmentを複数べん毛の主診断とし，close pair・participationは補助指標とする．旧Phase 2.7代表条件は後続のbasal freedom・dataset modelで置換されており，canonical model条件としては使用しない．
- **Evidence:** Task P2-7-006，`docs/phase2/phase2_7_flag_helix_axis_diagnostics.md`，`tests/test_helix_axis.py`，`tests/test_phase2_7_bundling_sweep.py`．

### P2-D07: 強いattach-frame固定ではなく，弱いposition-only補強を採用する

- **Status:** adopted
- **Background:** 後方条件ではhook過伸長が生じたためattach-frame補強を導入したが，強いposition・tangent固定はroot/body/flagellumの相対自由度を抑え，剛体回転様の挙動やproximal flag bond failureを生じさせた．
- **Change:** attach-first，body-axis angle，first-second，attach-frame position，attach-frame tangentを分離し，torque distributionとbody roll／axis-centered spinを比較した．
- **Comparison:** 強い`position=3`, `tangent=1.5`，first-second補強，torque distribution 2×2，position-only sweepを比較した．
- **Result:** strong tangent固定はroot軸まわりの相対spinも抑えやすく，2×2 torque distribution候補は長時間後方条件を解決しなかった．lateral条件では`no_frame`が自然な軸中心回転の参照となり，posterior条件では弱いposition-only補強`1.25`がhook安定性と相対回転の最良の折衷だった．
- **Interpretation:** 問題の主因はtorque profile単独ではなく，basal/root自由度を過度に固定することだった．
- **Decision:** `motor.local_attach_frame_position_scale=1.25`を2010 project defaultとする．tangent scaleとfirst-second spring scaleは`1.0`を維持し，torque distribution profileは`diffusive`を維持する．strong attach-frame候補とuniform／axis projection候補はcanonicalとして採用しない．
- **Evidence:** Tasks P2-8-082，P2-8-094，P2-8-097，P2-8-103，Issues #82，#94，#97，#103，ADR 0007．

---

## Analysis dataset・execution contract

### P2-D08: べん毛本数分析ではfeature registryとsample単位を固定する

- **Status:** adopted
- **Background:** RUN固定simulationから，べん毛本数差を再現可能に解析するには，feature名，QC，NaN，sample ID，timeseriesの契約が必要だった．
- **Change:** featureを`metadata`, `quality`, `cell_translation`, `cell_orientation`, `flagella_axis`, `cell_flagella_relation`, `diagnostics`へ分類し，YAML registryを正本とした．
- **Result:** `dataset_id`と`sample_id`を持つsummary・QC・sample別timeseriesを生成できるようになった．定義不能値は`NaN`とし，除外数を記録する契約を設けた．
- **Interpretation:** quality・diagnosticsはdatasetの妥当性確認に必要だが，原則としてML入力候補へ直接含めるべきではない．
- **Decision:** feature registry YAMLをfeature名・カテゴリの正本とする．quality・diagnosticsは原則ML入力候補から除き，`NaN`は集計時に欠損数を記録する．datasetはsummaryとsample別timeseriesを保持する．
- **Evidence:** Tasks P2-8-008，P2-8-009，P2-8-010，P2-8-014，`docs/phase2/phase2_8_flagella_count_feature_definitions.md`，`conf/phase2_analysis/flagella_count_behavior_features.yaml`．

### P2-D09: 付着点seedとhelix phase seedを分離する

- **Status:** adopted
- **Background:** 単一global seedだけでは，付着位置と初期helix phaseのばらつきを独立に制御・評価できなかった．
- **Change:** `seed.attach_seed`と`seed.phase_seed`を分離し，未指定時は`seed.global_seed`へfallbackする契約を追加した．`seeded_surface`とcenter-priority配置を導入した．
- **Comparison:** 同一seed再現性，attach seedのみ変更，phase seedのみ変更，legacy seed互換を比較した．
- **Result:** 付着点とhelix phaseを独立に再現できるようになった．小さいattach seedではbody中央層を優先し，その範囲を越えると制約なし配置へ進む．
- **Interpretation:** center-priorityは特定位置へ固定するものではなく，初期datasetで極端に偏った付着配置を避けるためのsampling policyである．
- **Decision:** attach seedとphase seedを独立軸とする．legacy seedは両seedへ同じ値を設定する互換形式として扱い，center-priorityの正確な展開数はbuilder・config・testを正本とする．
- **Evidence:** Tasks P2-8-012，P2-8-084b，P2-8-084c，`tests/test_model_builder.py`，`tests/test_params.py`．

### P2-D10: raw simulationとrender samplingを分離し，replay可能な出力を保持する

- **Status:** adopted
- **Background:** `dt_star=1.0e-4`では内部step数が多く，全stepを画像化すると実行時間とファイル数が大きくなる一方，保存時にstateを間引くと後続解析・再描画ができなくなる．
- **Change:** user-facingな複数条件実行を`run_multi_run.py`，診断gridを`run_sweep.py`，定量plotを`plot_heatmap.py`，保存済みstateの再描画をreplay CLIへ分離した．
- **Result:** `summary.csv`, `trajectory.csv`, `state_archive.npz`, `run_manifest.json`をreplayable output contractとして固定した．raw sampleは全step情報を保持し，軽量化はrender側のfps samplingで行う方針となった．保存段階のstrideは廃止し，canonical sweep名を`shape_stability_grid`，sweep summary名を`summary.csv`とした．通常診断は軽量な`run_summary.json`を先に読み，raw `step_summary.csv` は限定抽出時のみ使う．
- **Interpretation:** simulation情報を保存段階で不可逆に削減せず，描画・レビュー時にsamplingする方が，再現性と実行効率を両立できる．
- **Decision:** raw outputを解析・replay可能な形で保持し，通常renderでは全内部stepを描画しない．full-step renderは診断時のみ明示指定し，Issue番号付き使い捨てscriptをcanonical entrypointにしない．
- **Evidence:** Tasks P2-9-009，P2-8-078，P2-8-081，P2-8-084，P2-8-100，P2-8-DTSTAR，P2-SCRIPT-096，`conf/phase2_multi_run/README.md`，`scripts/README.md`．

### P2-D11: dataset v0はtraining datasetではなくdiagnostic baselineとする

- **Status:** diagnostic
- **Background:** 現行モデルでべん毛本数差が特徴量へ現れるかを確認するため，RUN固定の初期datasetが必要だった．
- **Change:** `n_flagella=[1,2,3,6]`，attach seed 3条件，phase seed 3条件の36 sampleを，motor torque `2.0e-20 N m`，1.0 sで生成した．generic multi-runとdataset builderを同じprofileへ統合した．
- **Result:** `n=1,2,3`は各9 sampleがstrict passし，speed，straightness，angular velocity，hook driftなどに本数差が見られた．`n=6`は全9 sampleでflag failureとなり，seed固定`n=4,5,6`でもflag failureとbody spring failureが併発した．
- **Interpretation:** v0はfeature候補とfailure modeを探索する価値はあるが，Phase 3/4 training datasetとして使用できる品質ではない．
- **Decision:** v0をdiagnostic baselineとして保持する．v0では`n=1,2,3`のみ分析対象とし，`n>=4`をtraining candidateに含めない．
- **Evidence:** Tasks P2-8-013，P2-8-015，P2-8-016，Issues #71，#113，#117，`conf/phase2_multi_run/flagella_count_behavior_v0.yaml`．

### P2-D12: dataset v1のtraining candidateはRUN固定`n=1,2,3`とする

- **Status:** adopted
- **Background:** v0では`n>=4`のflag・body failureが支配的だったため，flag springとbody stiffnessのproject-specific補強を比較し，改善モデルでdatasetを再生成する必要があった．
- **Change:** `stiffness_scales.flag_spring`と`stiffness_scales.body`を探索し，v1では`flag_spring=2.25`, `body=2.5`を使用した．全step gateをdataset QCへ反映した．
- **Comparison:** baseline `n=4,5,6`，flag spring・body scale候補，narrow sweep，v0とv1を比較した．
- **Result:** v1の`n=1,2,3`は27/27 sampleが全時間strict passだった．`n=4`は6/9 strict pass，3/9 flag failで，`n=5,6`は安定training候補を得られなかった．ユーザー定性評価で`n=1,2,3`が承認された．
- **Interpretation:** stiffness補強はn=4を改善したが，本数範囲全体の物理安定性を解決していない．途中failure後に最終stateで回復してもtraining candidateへ戻すべきではない．
- **Decision:** dataset v1のPhase 3/4 training candidateをRUN固定`n=1,2,3`とする．`n=4`はdiagnostic-only，`n>=5`はcanonical training scope外とし，training candidateには全時間strict passを要求する．
- **Evidence:** Tasks P2-8-018，P2-8-019，P2-8-020，Issues #115，#116，#119，`conf/phase2_multi_run/flagella_count_behavior_v1.yaml`．

### P2-D13: dataset versionとrevisionを分離し，configを条件の正本とする

- **Status:** adopted
- **Background:** config filename，dataset ID，output path，model条件が対応していないと，異なる物理モデルのdatasetを追跡できない．
- **Change:** model条件が変わるdataset versionと，同一model内の再生成・品質修正を表すrevisionを分離した．
- **Result:** canonical dataset IDを`v0`, `v1`のように管理し，同一version内の更新は`r1`, `r2`で表す．`v0.1`や`v0-1`は使用せず，詳細なmodel条件はconfigのmetadataとbase overridesを正本とする．旧dataset ID・config名はhistorical aliasとして記録する．
- **Interpretation:** dataset IDだけへ全条件を埋め込まず，version registryとconfig provenanceを組み合わせる方が追跡しやすい．
- **Decision:** dataset version・revision・model provenanceをregistry contractとして維持する．
- **Evidence:** Task P2-8-017，Issue #118，`docs/phase2/phase2_8_dataset_version_registry.md`，`conf/phase2_multi_run/README.md`．

### P2-D14: Phase 3へはdataset v1のRUN固定`n=1,2,3`を渡す

- **Status:** adopted
- **Background:** Phase 3のclip生成・tracking・metadata設計を開始するには，Phase 2から渡すbaseline，対象外条件，ID責務を固定する必要があった．
- **Change:** run，source video，track，clip，frameの単位と，各IDの責務を定義した．実動画と擬似動画を別入力経路から共通clip schemaへ収束させる方針とした．
- **Comparison:** 3D feature差がbody-centered XY投影後にも残るか，27独立runをseed group単位で評価した．
- **Result:** handoff baselineはdataset v1のRUN固定`n=1,2,3`とした．XY投影のgrouped nearest-centroid baselineは19/27，accuracy 0.704だった．body-centered renderではcell translationが見えないため，向き変化などを主要2D signalとして扱う必要がある．
- **Interpretation:** 2D上にも本数差の初期signalは残るが，0.704はtraining性能の保証ではない．clip長，独立run数，画像由来特徴はPhase 3/4で評価する．
- **Decision:** Phase 3 MVPへdataset v1の`n=1,2,3`を渡す．real videoはdetection/tracking経路，pseudo videoはGT passthrough経路とし，2D separability結果はfeasibility evidenceであってML acceptance gateにはしない．
- **Evidence:** Tasks P2-8-021，P2-8-022，Issues #125，#126，`docs/phase2/phase2_8_phase2_to_phase3_handoff.md`．

---

## Model correspondence・profiles

### P2-D15: 2010 projectのspringは`fene_fraenkel`，初期螺旋軸は後方をdefaultとする

- **Status:** adopted
- **Background:** 2010実装のspring式には，論文対応formulationと過去互換formulationが混在していた．また，project通常条件では初期螺旋軸の向きを明示する必要があった．
- **Change:** 論文対応の`fene_fraenkel`と過去互換`legacy`をselectorで分離し，motor-off `0.1 tau`とmotor-on RUN `1 tau`で比較した．2010 project profileの初期螺旋軸を菌体後方`0 deg`とした．
- **Result:** 両spring formulationは比較runを完走した．論文一致を優先して`fene_fraenkel`をdefaultとした．後方初期軸条件は0.5 sでfiniteかつshape gateを通過し，ユーザー定性評価でも承認された．
- **Interpretation:** `legacy`は過去互換性のため残すが，新しいcanonical runは論文対応formulationを用いる．後方初期軸はpaper profileではなくproject profileの通常条件である．
- **Decision:** 2010 projectのspring defaultを`fene_fraenkel`，`initial_helix_axis_from_rear_deg=0`とする．selector欠落時のlegacy compatibilityは維持し，2010 paper profileではpaper conditionを優先する．
- **Evidence:** Tasks P2-MODEL-163，P2-MODEL-171，Issues #163，#171，`docs/phase2/phase2_163_2010_potential_correspondence.md`，ADRs 0009，0011．

### P2-D16: 2010／2015のproject・paper profileとcanonical時間schemaを分離する

- **Status:** adopted
- **Background:** 2010 projectの互換条件，2010 paper条件，2015 project候補，2015 paper条件を1つのconfigで扱うと，モデル由来・解像度・実装状態・時間換算が不明確になる．
- **Change:** 4 profileを別configにし，`model_profile`とsource config pathをmanifestへ保存した．時間指定を`time.duration.value`, `time.duration.unit`, `time.integration.dt_star`へ統一した．
- **Result:** 2010 projectは15 body + 11×3 flagellum，2010 paperは15 body + 15×3 flagellum，2015 project/paperは30 body + 30×3 flagellumとした．2010 projectは`tau_s=1.0`互換を維持し，2010 paper／2015はreference torqueから`tau_s`を導出する．durationは`tau`と`s`を同一schemaで扱える．
- **Interpretation:** profileごとにgeometry，time scale，motor model，implementation statusを追跡する必要がある．
- **Decision:** default profileを`conf/sim_swim_2010.yaml`とする．legacy time keyはdeprecated compatibilityとして残し，新旧time keyの矛盾はerrorとする．canonical profileでは`motor.torque_Nm=-1` sentinelを使用せず，profile provenanceをmanifestへ保存する．
- **Evidence:** Tasks P2-MODEL-164，P2-MODEL-165，Issues #164，#165，ADRs 0010，0012．

### P2-D17: 2015 refined profileでは120-bead geometryとpaper-inspired body reactionを使用する

- **Status:** adopted
- **Background:** Kong et al. (2015)条件を評価するため，2010 triangular geometryとは別に，30 body beadsと30-bead flagellaを持つrefined geometryおよびbody reaction modelが必要だった．
- **Change:** 30 body beadsを正六角形×5層として構成し，flagellum 30 beads×3本と0.25 b Hookを実装した．2015 profileではbody diagonal braceをOFFとし，paper profile用に`hook_coupled_body_reaction`を実装した．
- **Comparison:** paper本文のpentagon記述と，Fig. 1・30 body beads・5 layersとの整合，project比較motorとpaper-inspired motorを比較した．
- **Result:** 六角形×5層が30 beadsおよび図と整合するproject inferenceとして採用された．paper motorはflagellum基部3 beadsへzero-net-force torque driveを与え，body側へ局所reactionを返す．project profileは`root_torque_segment_couples`，paper profileは`hook_coupled_body_reaction`とした．
- **Interpretation:** geometryとmotor reactionはいずれも論文の完全再現ではない箇所を含むため，paper値，figure inference，implementation assumption，project comparisonをprovenance上分ける必要がある．
- **Decision:** 2015 refined geometryを採用し，六角形×5層をproject inference，paper motorを`paper_inspired_approximation`として記録する．2015 projectは後方初期軸，paper profileはpaper conditionを維持し，supported昇格はStage A評価後に判断する．
- **Evidence:** Tasks P2-MODEL-166，P2-MODEL-167，Issues #166，#167，`docs/phase2/phase2_167_2015_paper_conditions.md`，ADRs 0013，0014．

---

## Pending decisions

### P2-D18: 2015 refined modelのStage A採否を確定する

- **Status:** pending
- **Background:** 2015 refined geometryとmotor dynamicsは実装済みだが，短時間安定性，integration step感度，定性挙動を確認しない限りsupported profileへ昇格できない．
- **Change:** motor-off／motor-on，project／paper，`dt_star=1e-4`／`1e-5`を比較できるStage A runnerと判定出力を用意した．
- **Current result:** motor-off `0.1 tau`とcanonical motor-on `1 tau`は完走した。`dt_star=1e-4`はproject profileの先頭`0.1 tau`ペアでは回転・姿勢・bead差の基準を満たしたが，paper profileでは回転量差が10%を超えた。project後方束化のtorque scale `0.5, 1, 2` × `dt_star=1e-5, 1e-4`は全6 conditionが完走し、有限性・shape・motor action-reaction gateを通過した。grid replayのユーザー目視では明瞭なcollapse / fly-awayは認められなかった。
- **Interpretation:** Stage Aの短時間検証と定性可視化は完了したが、`dt_star=1e-4`をproject/paper共通defaultへ昇格させる根拠はない。短い実時間では小さな誤差を過大評価し得るため、同一実時間での定量評価と計算効率評価を分けて行う。
- **Decision:** canonical `dt_star=1e-5`を維持し，`dt_star=1e-4`は非canonical referenceとする。本PRのStage A検証は完了とする。`dt_star`の有効性説明はIssue #61、reference torque policyはIssue #183、dataset v2採択向けのtorque・`dt_star`・べん毛数検証はIssue #184で実施し、2015 profileのsupported採否はそれらの後に判断する。
- **Evidence:** Issue #168，PR #176，Issues #61・#183・#184，`docs/phase2/phase2_168_2015_stage_a_validation.md`，parent Issue #154．

### P2-D19: 大量step runはcompact output policyを明示選択する

- **Decision:** 既存短時間 workflow は `output.policy: debug` の全step state/CSV 互換を維持する。長時間候補は `compact` を明示し、全内部step QC をオンライン集約しながら state archive を物理時間一様に保存する。標準 cadence は 1 ms とする。compact archive は replay/通常解析用であり、未定義の全step将来指標の完全再構成を保証しない。
- **Evidence:** Issue #186，`conf/phase2_multi_run/2015_project_compact_0p5s.yaml`，`docs/phase2/phase2_run_summary_contract.md`．

### P2-D19: dataset v2前にn=3長時間非定常回転とproximal failureを診断する

- **Status:** pending
- **Background:** dataset v1 r1の3.0 s RUN固定`n_flagella=3`で，非定常回転とproximal flag bond failureが確認された．原因未確定のままdataset v2を生成すると，同じfailureを長時間datasetへ持ち込む可能性がある．
- **Change:** 既存raw outputから，first-fail前後のangular velocity，speed，axis alignment，proximal bond，contact／penetration proxyを同期解析する計画とした．
- **Current result:** 診断基盤・原因仮説・解決策候補の整理が進行中であり，物理モデル変更は未採択である．
- **Interpretation:** dataset v2 model変更は，failure位置・seed依存・contact proxyとの時間関係を確認してから判断すべきである．
- **Decision:** Issue #158完了までは原因を断定しない．dataset v2 coreはRUN固定`n=1,2,3`を対象とし，training candidateには全時間strict passを要求する．source durationは5 sを候補とし，attach seed×phase seedを独立runとして扱う．
- **Evidence:** Task P2-8-023，Issues #157，#158，`docs/phase2/phase2_158_v1_r1_nf3_proximal_diagnostics.md`．

### P2-D20: RUN–TUMBLEはRUN dataset core完了後に段階実装する

- **Status:** pending
- **Background:** 現行canonical datasetは`motor.enable_switching=false`のRUN固定であり，motor reversalやpolymorph switchingを含まない．
- **Change:** Tumble導入を，設計，motor reversal，polymorph switching，run-and-tumble評価へ分割した．
- **Interpretation:** motor reversalとnormal／semicoiled／curly1切替を一度に導入すると，failure原因と受入条件を分離できない．RUN baselineを壊さない段階実装が必要である．
- **Decision:** dataset v2 RUN core完了後にIssue #69で扱う．各stageを別review resultで評価し，RUN固定挙動の回帰を必須とする．state transitionはstep summaryまたは専用CSVへ記録する．
- **Evidence:** Task P2-9-010，Issue #69．
