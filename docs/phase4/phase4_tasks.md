# Phase 4 Decisions

この文書は，Phase 4の後続開発で必要となる採択判断と根拠を保持する．

Issue単位の進捗台帳ではない．現在の作業状態は`phase4_current.md`，datasetの条件と実行方法はconfigおよび`scripts/README.md`を正本とする．

## P4-D01: Phase 3 common clip datasetのloader contract

* **Status:** adopted
* **Background:** Phase 4がPhase 3出力を再検出せず利用するため，clip，metadata，splitの整合性を検査する必要があった．
* **Change:** `manifest.json`，`clip_metadata.jsonl`，`split_summary.csv`，`.npy` clipを読むloaderとaudit helperを実装した．
* **Result:** clipのdtype・shape・frame数，label，split，`group_key` leakageを検出できるようになった．
* **Decision:** Phase 4はPhase 3 common clip datasetを直接読み込み，独自の再検出や再分割を行わない．
* **Evidence:** Issue #146，PR #147，`src/flagella_estimation/phase4/dataset.py`，`tests/test_phase4_clip_dataset_loader.py`

## P4-D02: 診断baseline classifier

* **Status:** diagnostic
* **Background:** 学習・評価pipelineとdataset contractを確認するため，軽量で決定的なbaselineが必要だった．
* **Change:** clipの輝度，foreground，時間差分，重心移動，radial spreadを特徴量とし，standardization付きnearest-centroid classifierを実装した．
* **Result:** Phase 3のgrouped splitを維持したまま，学習，予測，metrics，confusion matrix，model artifactを生成できた．
* **Interpretation:** 少数のpseudo groupsに対する性能は，最終modelの性能や実動画への一般化を示さない．
* **Decision:** このbaselineはpipeline診断とdataset比較に使用し，最終model候補として採択しない．
* **Evidence:** Issue #148，PR #149，`conf/phase4/baseline_v1.yaml`，`tests/test_phase4_baseline_classifier.py`

## P4-D03: grouped learning curve

* **Status:** diagnostic
* **Background:** 同一runから複数clipを生成しても独立情報は増えないため，clip数ではなく独立group数に基づくlearning curveが必要だった．
* **Change:** 同一`group_key`のclip特徴を平均し，1 groupを1 feature vectorとして評価するevaluatorを実装した．
* **Comparison:** train+validation pool内でclass-balancedなtrain / holdout groupsを作り，test splitを保護した．
* **Result:** 現行pseudo-v1と診断baselineでは，4 training groups/class以降でmacro F1のplateauを確認した．
* **Interpretation:** protected testが1 group/classであり，実動画domainも未評価であるため，一般的な必要run数を4とは決められない．
* **Decision:** `k=4`はpseudo-v1内の診断結果として扱い，必要run数の採択範囲とprotected evaluation条件は#129で決める．
* **Evidence:** Issue #150，PR #151，Issue #129，`conf/phase4/grouped_learning_curve_v1.yaml`，`tests/test_phase4_learning_curve.py`

## P4-D04: dataset freezeとmixing policy

* **Status:** adopted
* **Background:** dataset version，物理model，render条件，behavior条件の異なるdataが無条件に混在すると，dataset解釈と評価が不明確になる．
* **Change:** dataset，model，render，class，clip，QC，group，Phase 2 source provenanceを検査するmachine-readable freeze auditを実装した．
* **Result:** trainingとgrouped learning curveが同じstrict validatorを使用し，source manifestと各runの解決済みconfigまで追跡できるようになった．
* **Decision:** dataset v1 baselineはRUN固定，`n_flagella=1,2,3`，baseline torque，Brownianなしとする．条件混在とversion境界はADR 0015を正本とする．
* **Evidence:** Issue #128，PR #152，`docs/adr/0015_dataset_mixing_and_versioning.md`，`conf/phase4/dataset_freeze_v1.yaml`，`conf/phase4/dataset_freeze_v1_r1.yaml`，`tests/test_phase4_dataset_freeze.py`

## P4-D05: 3秒duration / seed study

* **Status:** adopted
* **Background:** 1秒sourceだけでは，clip開始時刻，初期過渡状態，attach / phase seed差を十分に比較できなかった．
* **Change:** dataset v1と同じ物理条件を3秒へ延長し，`n_flagella=1,2,3`，attach seed 3条件，phase seed 3条件の27 independent runsを評価した．
* **Comparison:** 0.25秒，0.5秒，1.0秒windowの3D / 2D特徴，run / window QC，attach×phase seed差を比較した．
* **Result:** 0.5秒と1.0秒で主要2D特徴のクラス順位は概ね維持された．strict QCは`n=1: 9/9`，`n=2: 9/9`，`n=3: 5/9`であり，`n=3`はattach seed依存が大きかった．
* **Interpretation:** 0.5秒clipは現行baselineとして維持できるが，`n=3`の4 fail runを正式なv2 training dataへ含めることはできない．
* **Decision:** MVP defaultを0.5秒のまま維持する．v1 r1ではwindow QCとpre-first-fail provenanceを保持する。#158の診断は現行v1 r2 campaignのphysical-failure blockerではなく，非定常回転とseed差はPhase 3でwithin-class variationとして評価する。今後のdataset sourceとfreezeは#157およびPhase 3 Gateへ分離する．
* **Evidence:** Issue #129，Issues #155・#157・#158，`conf/phase2_multi_run/flagella_count_duration_3s_r1.yaml`，`scripts/04_phase4/run_duration_seed_study.sh`，`tests/test_phase4_duration_study.py`

## P4-D06: warmupとearly clip

* **Status:** pending
* **Background:** simulation開始直後のmechanical relaxationが，class分離性と学習結果へ影響する可能性がある．
* **Change:** early clipを削除せず，全artifactと時刻情報をPhase 4へ渡せるようにした．
* **Comparison:** `warmup_s=0 / 0.5 / 1.0`について，時間帯別特徴，grouped performance，run内・run間分散を比較する．
* **Decision:** early clipをtraining candidate，augmentation，diagnostic-onlyのどれとして扱うかは#155で決める．
* **Evidence:** Issue #155，Issue #159，`conf/phase4/dataset_freeze_v1_r1.yaml`

## P4-D07: dataset v2とRUN-TUMBLE

* **Status:** pending
* **Background:** dataset v1 r1には`n_flagella=3`の長時間failureがあり，RUN-TUMBLEも現行RUN固定baselineとは異なるbehavior regimeである．v1 r1 failure診断は現行v1 r2 campaignのphysical-failure blockerではなく，非定常回転とseed差はPhase 3のwithin-class variationとして評価する．
* **Change:** dataset v2をRUN固定coreとRUN-TUMBLE scopeへ分離して進める方針を定めた．
* **Decision:** #157でcanonical sourceのpilotを生成し，Phase 3 Gate後にRUN固定coreをfreezeする．RUN-TUMBLEは#69とv2 core通過後に#145で別`behavior_profile` / `dataset_scope`として追加する．
* **Evidence:** Issues #145，#157，#158，ADR 0015
