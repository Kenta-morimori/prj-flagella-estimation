# Phase 3 Current

この文書は，Phase 3作業の入口となる短い現在地ドキュメントである．

## Goal

Phase 3の目的は，実顕微鏡動画とPhase 2擬似動画を，Phase 4で共通利用できる個体clipとmetadataへ変換することである．

## Current Baseline

* 実動画と擬似動画で共通のclip metadata schemaを採用している．
* Phase 2擬似dataはGround Truth passthroughにより，再検出せずclip datasetへ変換できる．
* dataset v1 r1の3秒runから，0.5秒non-overlap clipを生成するpipelineを実装済みである．
* Issue #199のv1 r2 follow-on条件では，2秒sourceを`0.25/0.5/1.0/1.5 s` non-overlap windowへ分割し，n=4をdiagnostic-onlyとして保持する準備を完了した．
* canonical renderは`body_capsule_orthographic_v1`である．
* Phase 4のloader，baseline，grouped learning curve，freeze auditへ接続済みである．
* 実動画のdetection / tracking経路は未実装である．

## Active Work

* #199: tau-linked `T=2.5e-20 N m`, `dt*=1e-3`, body scale 1.0でv1 r2 sourceをユーザー実行し，3D+2D feature comparisonで本数差を確認する．

## Next Queue

1. #8: 実顕微鏡動画の入力条件と必要metadataを整理する．
2. #9: 実菌体サイズとscale normalizationの基準を整理する．
3. #17: 実動画条件に基づくrender profileを定義する．

## Blockers

* 実動画のdetection / tracking実装は，#8と#9の入力・scale条件に依存する．
* 実動画domainに合わせたrender variationは，#8の撮影条件確定まで採択しない．

## Context Routing

* 採択判断: `docs/phase3/phase3_tasks.md`
* Human-readable schema: `docs/phase3/phase3_1_clip_metadata_schema.md`
* Machine-readable schema: `schemas/phase3_clip_metadata.schema.json`
* Dataset mixing policy: `docs/adr/0015_dataset_mixing_and_versioning.md`
* CLIと実行方法: `scripts/README.md`
