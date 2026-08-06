# Phase 4 Current

この文書は，Phase 4作業の入口となる短い現在地ドキュメントである．

## Goal

Phase 4の目的は，Phase 3 common clip datasetを読み込み，べん毛数`n_flagella`の学習・評価を行うことである．

## Current Baseline

* Phase 3 common clip datasetのloader，診断baseline，grouped learning curve，dataset freeze auditを実装済みである．
* dataset v1 r1の0.5秒non-overlap clipを現行baselineとする．
* training candidateはRUN固定の`n_flagella=1,2,3`である．
* 同一run由来のclipは`group_key`で束ね，独立sample数へ重複計上しない．
* baseline classifierの性能はpipeline診断値であり，最終model採択や実動画への一般化を意味しない．
* 3秒full-factorial 27 runの解析とPhase 3 clip datasetへの変換は完了している．
* 3秒runのstrict QCは`n=1: 9/9`，`n=2: 9/9`，`n=3: 5/9`であり，`n=3`の物理品質改善が必要である．

## Active Work

* #129: pseudo-v1に必要な独立run数の採用範囲とprotected evaluation条件を決める．

現行pseudo-v1と診断baselineでは，4 training groups/classからlearning curveがplateauした．ただし，protected testは1 group/classであり，一般的な必要run数を4と結論づける根拠にはしない．

## Next Queue

1. #155: `warmup_s=0 / 0.5 / 1.0`を比較し，early clipの用途を決める．
2. #158: `n_flagella=3`の非定常回転とproximal failureを診断する．
3. #157: 修正済み物理modelによるdataset v2 RUN固定coreを生成・freezeする．
4. #145: v2 core通過後にRUN-TUMBLE scopeを追加する．

## Blockers

* 必要run数の一般化には，追加のprotected evaluation groupsが必要である．
* dataset v2のfreezeには，`n_flagella=3`の長時間failureの説明と解消が必要である．
* RUN-TUMBLE datasetは，#69とdataset v2 RUN固定coreの完了に依存する．
* 実動画への一般化評価は，Phase 3の実動画detection / tracking経路に依存する．

## Context Routing

* 採択判断: `docs/phase4/phase4_tasks.md`
* Phase 3の現在地: `docs/phase3/phase3_current.md`
* Common clip contract: `docs/phase3/phase3_1_clip_metadata_schema.md`
* Dataset mixing policy: `docs/adr/0015_dataset_mixing_and_versioning.md`
* Dataset freeze config: `conf/phase4/dataset_freeze_v1_r1.yaml`
* CLIと実行方法: `scripts/README.md`
