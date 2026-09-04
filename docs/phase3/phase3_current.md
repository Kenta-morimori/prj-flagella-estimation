# Phase 3 Current

この文書は，Phase 3作業の入口となる短い現在地ドキュメントである．

## Goal

Phase 3の目的は，canonical physical simulation modelからべん毛本数による3D運動差と2D観測可能性を評価し，実顕微鏡条件に対応した個体clip / metadata datasetをfreezeしてPhase 4へ渡すことである．
流れは，controlled study，3D特徴量，ideal 2D projection，pseudo microscopy，observable feature，robustness，dataset freezeとする．

## Current Baseline

* 実動画と擬似動画で共通のclip metadata schemaを採用している．
* Phase 2擬似dataはGround Truth passthroughにより，再検出せずclip datasetへ変換できる．
* dataset v1 r1の3秒runから，0.5秒non-overlap clipを生成するpipelineを実装済みである．
* Issue #199のv1 r2 follow-on条件では，2秒sourceを`0.25/0.5/1.0/1.5 s` non-overlap windowへ分割し，n=4をdiagnostic-onlyとして保持する準備を完了した．
* canonical renderは`body_capsule_orthographic_v1`である．
* Phase 4のloader，baseline，grouped learning curve，freeze auditへ接続済みである．
* 実動画のdetection / tracking経路は未実装である．
* v1 / v1 r1 / v1 r2のseed sweepは探索的baselineとして保持し，canonical Phase 3 GateはPhase 2 canonical model freeze後に再評価する．

## Active Work

* #205のcanonical physical model freezeに先行して，#204のfeature contractと#8 / #9の実顕微鏡条件を整理する．

## Next Queue

1. #204 / #126: 3D・2D特徴量を固定し，controlled count effectとideal 2D observabilityを評価する．
2. #8 / #9 / #17: 実動画条件と擬似顕微鏡forward modelを定義する．
3. #157: canonical sourceのpilot，Phase 3 Gate後のdataset v2 freeze，Phase 4 handoffを行う．

## Blockers

* 実動画のdetection / tracking実装は，#8と#9の入力・scale条件に依存する．
* pseudo microscopyのrealistic profileとobservation robustnessは，#8 / #9の根拠確定までfreezeしない．
* dataset v2 final freezeは，controlled effect，ideal / realistic 2D observability，限定robustnessのPhase 3 Gateに依存する．

## Context Routing

* 採択判断: `docs/phase3/phase3_tasks.md`
* Human-readable schema: `docs/phase3/phase3_1_clip_metadata_schema.md`
* Machine-readable schema: `schemas/phase3_clip_metadata.schema.json`
* Dataset mixing policy: `docs/adr/0015_dataset_mixing_and_versioning.md`
* CLIと実行方法: `scripts/README.md`
