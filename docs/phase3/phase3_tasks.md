# Phase 3 Decisions

この文書は，Phase 3の後続開発で必要となる採択判断と根拠を保持する．

Issue単位の進捗台帳ではない．現在の作業状態は`phase3_current.md`，詳細なschemaは対応するschema文書とmachine-readable schemaを正本とする．

## P3-D01: 実動画・擬似動画の共通clip schema

* **Status:** adopted
* **Background:** 実動画はdetection / trackingを必要とする一方，Phase 2擬似動画はGround Truth trackを利用できる．入力経路が異なっても，Phase 4へ渡す形式は共通化する必要があった．
* **Change:** `detection_tracking`と`gt_passthrough`を同じclip metadata schemaへ収束させた．経路差は`processing_mode`，`label_source`，provenanceで表現する．
* **Result:** machine-readable schema，最小fixture，contract testを追加し，必須fieldと条件付きvalidationを固定した．
* **Decision:** Phase 4は入力経路ではなく共通schemaを読み，dataset splitには`track.group_key`を使用する．
* **Evidence:** Issue #127，PR #142，`docs/phase3/phase3_1_clip_metadata_schema.md`，`schemas/phase3_clip_metadata.schema.json`，`tests/test_phase3_clip_metadata_schema.py`

## P3-D02: clip時間長と独立sample単位

* **Status:** adopted
* **Background:** 同一runやtrackから複数clipを生成できるため，clip数と独立sample数を分離する必要があった．
* **Change:** `0.25 s`，`0.5 s`，`1.0 s`とnon-overlap / overlapを比較し，split単位を`track.group_key`へ統一した．
* **Comparison:** `0.25 s`は短時間耐性，`1.0 s`はfull-run比較，overlapは反復観測として扱った．
* **Result:** `0.5 s`で現行sourceから複数windowを生成でき，Phase 4の入力単位として利用できた．
* **Interpretation:** 同一run内clip，overlap window，同一raw runのrender variationは独立sampleではない．
* **Decision:** MVP defaultを`0.5 s`，window policyを`non_overlap`とする．独立数はpseudo dataでは`run_id`，実動画では`source_video_id`と`track_id`の組み合わせで数える．
* **Evidence:** Issue #129，`conf/phase3/gt_passthrough_v1.yaml`，`conf/phase3/gt_passthrough_v1_r1_duration_3s_clips.yaml`

## P3-D03: Phase 2擬似dataのGT passthrough

* **Status:** adopted
* **Background:** Phase 2擬似dataではbody positionとlabelを取得できるため，再度detectionを行う必要がなかった．
* **Change:** state archiveからwindow生成，rasterize，metadata生成，grouped split，QC summaryまでを行うPhase 3 pipelineを実装した．
* **Result:** `.npy` clipと共通metadataを生成し，Phase 4が再検出せず読み込めるようになった．
* **Decision:** pseudo dataはGT passthroughを標準経路とし，将来の実動画経路はadapter境界の後段で同じmetadata builderへ接続する．
* **Evidence:** Issue #6，PR #144，`src/flagella_estimation/phase3/`，`scripts/03_phase3/build_clip_dataset.py`，`tests/test_phase3_gt_passthrough_pipeline.py`

## P3-D04: dataset v1 r1 common clip dataset

* **Status:** adopted
* **Background:** dataset v1 r1の3秒runを，Phase 4で扱える固定長clipへ変換する必要があった．
* **Change:** 27 independent runsを0.5秒non-overlap windowへ分割し，run / window QC，first-fail時刻，training candidate，diagnostic-onlyをmetadataへ伝搬した．
* **Comparison:** 3D / 2D replayにより，斜め姿勢とz方向へ正対した場合のsilhouetteを確認した．
* **Result:** 各runから5 clipを生成でき，`body_capsule_orthographic_v1`によるrigid capsule renderが目視評価で採択された．
* **Interpretation:** first-failを含むwindowとそれ以後は診断用であり，early clipは削除せずPhase 4側の選択条件として扱える．
* **Decision:** 全windowをartifactとして保持し，first-failを含むwindowとそれ以後をdiagnostic-onlyとする．canonical renderは`body_capsule_orthographic_v1`とする．
* **Evidence:** Issue #159，PR #161，`conf/phase3/gt_passthrough_v1_r1_duration_3s_clips.yaml`，Phase 3 replay，Phase 4 freeze audit

## P3-D05: dataset mixingとversion境界

* **Status:** adopted
* **Background:** render variation，seed違い，物理条件変更，model変更を同一datasetへ混在させると，split leakageとdataset解釈の不整合が生じる．
* **Change:** 条件変更を観測augmentation，独立run，domain variation，dataset version変更へ分類し，identifierとfreeze gateの責務を定めた．
* **Result:** Phase 4のmachine-readable freeze auditでdataset，model，render，group，source provenanceを検査できるようになった．
* **Decision:** 同一raw run由来のvariationは同じ`group_key`を維持する．domain variationやmodel変更をbaselineへ無条件に混ぜない．詳細はADR 0015を正本とする．
* **Evidence:** Issue #128，`docs/adr/0015_dataset_mixing_and_versioning.md`，`conf/phase4/dataset_freeze_v1.yaml`，`conf/phase4/dataset_freeze_v1_r1.yaml`

## P3-D06: 実動画のdetection / tracking経路

* **Status:** pending
* **Background:** 実動画にはPhase 2のGround Truth trackと`n_flagella` labelが存在せず，低コントラスト菌体と背景artifactを区別する必要がある．
* **Change:** 共通schemaにbody axis，長短径，detection confidence，tracking gapなどの任意fieldを定義した．
* **Result:** 初期のthresholdingとconnected componentsでは，小さい黒点，背景noise，リング状artifactを多く検出し，採用手法を確定できなかった．
* **Interpretation:** detectorとnormalizationを確定するには，実動画の入力条件と実菌体scaleの整理が先に必要である．
* **Decision:** #8で入力条件，#9でscale normalization，#17でrender条件を整理した後に実動画経路を実装する．共通出力schemaは変更しない．
* **Evidence:** Issues #6，#8，#9，#17，`docs/phase3/phase3_1_clip_metadata_schema.md`
