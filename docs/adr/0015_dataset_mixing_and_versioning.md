# ADR 0015: 学習datasetの条件混在とversion境界を明示する

* Status: Accepted
* Date: 2026-08-06
* Issues: #128, #129, #133, #145, #157

## Context

Phase 2由来datasetには，観測条件の変更，同じ物理条件からの独立run，物理domainの変更，modelやlabel範囲の変更が含まれる．

これらを同一datasetへ無条件に混ぜると，split leakage，provenance欠損，dataset解釈の不整合が生じる．

## Decision

条件変更を次の4種類へ分類する．

| class             | examples                                                    | rule                                                   |
| ----------------- | ----------------------------------------------------------- | ------------------------------------------------------ |
| 観測augmentation    | crop，rotation，brightness，noise，blur，codec，render profile    | label-preservingな範囲だけ使用し，同一raw sampleの`group_key`を維持する |
| 独立run             | attach seed，phase seed，initial pose，独立simulation seed       | 同一model・behavior regimeなら別groupとして追加できる                |
| domain variation  | torque，viscosity，Brownian，RUN-TUMBLE，motion evidenceを変えるfps | 根拠なしにbaseline trainingへ混ぜない                            |
| dataset version変更 | equation，constraint，stiffness，geometry，model，label range    | `dataset_version`と通常`model_id`を更新する                    |

## Identifier Responsibilities

* `model_id`: 物理modelと主要なmodel条件を識別する．
* `dataset_version`: datasetの論理解釈と採択範囲を識別する．
* `dataset_revision`: 同じversion内の再生成または条件追加を識別する．
* `render_id`: 同一raw runから生成した観測条件を識別する．
* `behavior_profile`: RUN固定やRUN-TUMBLEなどの挙動条件を識別する．
* `dataset_scope`: core，diagnostic，robustnessなどの用途を識別する．
* `group_key`: split leakage防止の正本とする．

同一raw runから生成したrender variation，overlap window，normalization variationは同じ`group_key`を維持する．

## MVP Dataset v1

MVP training baselineは次の条件とする．

* RUN固定
* `n_flagella=1,2,3`
* baseline torque
* Brownianなし
* `n_flagella=4`はdiagnostic-only
* 根拠のある小範囲の観測augmentationのみ許可

torque variation，Brownian，RUN-TUMBLE，model変更，label範囲変更をbaselineへ無条件に混ぜない．

training前にmachine-readable freeze auditを実行し，dataset，model，render，group，QC，source provenanceを検証する．

## Dataset v2

dataset v2では，修正済み物理modelによるRUN固定coreを先にfreezeする．

RUN-TUMBLEは別の`behavior_profile`または`dataset_scope`として追加し，v1 baselineやv2 RUN coreへ自動的に混ぜない．

評価は少なくとも次を分離する．

* v1 robustness
* v2 RUN core
* v2 RUN-TUMBLE
* 明示的に承認した混合training

## Consequences

* clip数と独立sample数を分離して扱う必要がある．
* 同一run由来variationはsplitをまたげない．
* domain variation追加時はfreeze policyとprovenanceを更新する必要がある．
* dataset version変更は過去datasetとの自動互換を意味しない．
* diagnostic-only条件は正式なtrain / validation / testへ含めない．

## Evidence

* Issue #128
* `conf/phase4/dataset_freeze_v1.yaml`
* `conf/phase4/dataset_freeze_v1_r1.yaml`
* `schemas/phase3_clip_metadata.schema.json`
* Phase 3 / 4 dataset freeze tests
