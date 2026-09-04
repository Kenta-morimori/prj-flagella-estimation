# ADR 0020: Phase 2–3の責務境界と研究Gate

## Status

Accepted

## Context

従来のPhase 2には3D物理simulation，2D擬似動画，dataset handoffが混在し，Phase 3はclip生成中心に記述されていた．また，canonical model selection Issue #205が3D / 2D識別性の評価を入力に含めており，freeze済みmodelを前提とするPhase 3評価と循環する余地があった．

## Decision

Phase 2は物理simulation modelの構築，論文・実装対応，数値安定性，物理妥当性，長時間安定性を担当し，その証拠だけでcanonical physical simulation modelをfreezeする．

Phase 3はfreeze済みmodelを入力として，controlled `n_flagella × attach_seed × phase_seed` study，3D feature，ideal 2D projection，realistic pseudo microscopy，pixel-observable feature，限定nuisance robustnessを評価する．Phase 3はdataset v2をfreezeしてPhase 4へhand offする．

Phase 3 Gateはcount effect，3D→ideal 2D information retention，realistic observation，body geometry / per-motor torque / viewing orientation / image quality bundleへのrobustness，provenanceを伴うdataset freezeとする．分類器の性能閾値とp値はGateにしない．MLによる未知groupへの予測性能はPhase 4の責務である．

既存v1 / v1 r1 / v1 r2のseed sweepは，Phase 2 QC履歴とPhase 3設計の探索的根拠として保持するが，canonical Phase 3 Gateの結果には使用しない．

## Consequences

- #205は#204/#126を待たずにPhase 2の選定Gateとして完結する．
- #126と#17はPhase 3の観測可能性評価を段階化する．
- #157はcanonical sourceのpilotをPhase 3評価の入力として先行生成できるが，final dataset freezeはPhase 3 Gateの後に行う．
- #225の流体診断は研究上有用でもPhase 3開始のhard blockerではない．
- v1 r1の`n_flagella=3` failure診断は現行v1 r2 campaignのphysical-failure gateではない。非定常回転とseed差は，Phase 3でwithin-class variationとして再評価する．
