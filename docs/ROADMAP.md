# Project Roadmap

このドキュメントは，プロジェクト全体の現在地，Issue階層，次に進める順序，ユーザー実行が必要な作業を一画面で追うための地図である。

詳細な実験履歴は `docs/phase2/phase2_current.md` / `docs/phase2/phase2_tasks.md` / `docs/codex-runs/*/review_result.json` を正本とする。このファイルでは，流れを把握するために主要Issueだけを束ねる。

## Reporting Unit

Codexが一気に実装を進める場合，報告単位は次の順にする。

1. PR単位: コード・docs・testsがまとまり，push済みの状態。
2. Review result単位: `docs/codex-runs/<run-id>/review_result.json` が `PASS` または診断価値のある `FAIL` になった状態。
3. User-run request単位: 長時間simulation，動画render，目視レビューなど，ユーザー実行または判断が必要になった状態。

報告には必ず以下を含める。

- こちらで完了したこと
- こちらで実行したコマンド
- こちらでは実行していないこと
- ユーザー実行が必要か
- ユーザー実行が必要な場合の exact command
- 確認すべき出力
- 判断ポイント
- 次に進める項目

原則として，ユーザー判断が不要な実装・軽量テスト・docs更新・PR作成まではCodex側で進める。判断が必要な部分だけを止める。

ユーザー判断で止める代表例:

- 長時間simulation / sweep / render の実行可否
- 2D/3D動画の自然さ，body deformation，helical shape preservation などの目視レビュー
- Phase 3/4 の仕様分岐，datasetに混ぜる条件，実データ取得条件
- PRをreadyにするか，mergeするか

## Current Top-Level Flow

```text
Project MVP
├─ Phase 2 RUN baseline and pseudo-video readiness
├─ Phase 3 common clip pipeline
├─ Phase 4 flagella-count training baseline
├─ Phase 2 physics extensions
└─ Tooling / workflow
```

## Current Active PRs

| PR | State | Purpose | Next |
| ---: | --- | --- | --- |
| #132 | open | プロジェクト roadmap と issue hierarchy 整理 | review gate / CI 通過後に merge |

## Issue Hierarchy

### Project / Phase Parents

| Issue | State | Role | Notes |
| ---: | --- | --- | --- |
| #10 | open | Phase 2 physics parent | 論文モデル実装・拡張の親Issue。物理モデル改善，Tumble，Brownian，安定化候補を束ねる。 |
| #71 | closed | Phase 2 RUN count-analysis parent | 3D RUN固定で `n_flagella` 差をdatasetとして評価した中Issue。#125でPhase 3へhandoff済み。 |
| #6 | open | Phase 3 clip pipeline parent | 実動画・擬似動画から共通個体clipを生成する実装親Issue。 |

### Phase 2 RUN Baseline / Dataset Handoff

| Issue | State | Parent | Role |
| ---: | --- | ---: | --- |
| #65 | closed | #71 | べん毛数分析用特徴量定義 |
| #72 | closed | #71 | 複数条件実行とdataset構築 |
| #73 | closed | #71 | 特徴量分布可視化 |
| #76 | closed | #10 / #71 context | attach_seed / phase_seed の分離 |
| #78 | closed | #71 | dataset replay render のfps指定 |
| #81 | closed | #71 | sweep実行時間短縮 |
| #112 | closed | #71 | 診断dataset v0条件 |
| #113 | closed | #71 | `n>=4` seed固定破綻境界 |
| #117 | closed | #71 | diagnostic dataset v0統計レポート |
| #118 | closed | #71 | `model_id` / dataset config命名 |
| #119 | closed | #71 | 改善モデル dataset v1生成 |
| #125 | closed | #71 | Phase 2→3 handoff baseline |
| #126 | open | #71 | 2D投影後の識別性確認 |

Current baseline:

- dataset: `outputs/phase2_analysis/flagella_count_behavior/datasets/v1`
- training candidate: RUN固定 `n_flagella=1,2,3`
- diagnostic-only: `n_flagella=4`
- out of scope for MVP: `n_flagella>=5`, TUMBLE, Brownian, torque variation, model changes

### Phase 3 / Phase 4 Handoff

| Issue | State | Parent | Role |
| ---: | --- | ---: | --- |
| #8 | open | #6 context | 実顕微鏡動画の取得・入力条件 |
| #9 | open | #6 context | 実顕微鏡像から菌体長さ特徴を分析 |
| #17 | open | Phase 2→3 context | 実顕微鏡像に合わせた擬似顕微鏡描画条件 |
| #127 | open | #6 | 実動画・擬似動画の共通clip / metadata schema |
| #128 | open | Phase 2→4 context | 学習datasetへ混ぜてよい条件変更 |
| #129 | open | Phase 3→4 context | 1 clip時間長と必要な独立run数 |

Recommended order:

1. #126: trajectory / frame / clip の2D識別性を確認する。
2. #127: common clip / metadata schema を固定する。
3. #129: clip時間長と独立run数を評価する。
4. #128: augmentation / domain variation / dataset version規則を固定する。
5. #6: Phase 3 clip pipeline を実装する。

### Phase 2 Physics Extensions

| Issue | State | Parent | Role |
| ---: | --- | ---: | --- |
| #82 | closed | #10 | hook過伸長対策の中Issue |
| #94 | closed | #82 | hook補強時のflag bond過伸長診断 |
| #97 | closed | #82 | torque distribution見直し |
| #100 | closed | #82 | sweep / multi-run / replay導線整理 |
| #103 | closed | #82 | basal freedom / attach-frame剛体回転診断 |
| #124 | open | #10 | `n>=4` 根元付近過伸長改善 |
| #69 | open | #10 | Tumble状態の段階実装 |
| #15 | open | #10 | Brownian項追加 |
| #41 | open | #10 | 参照論文 contour length 矛盾対応 |
| #61 | open | #10 context | `dt_star` 有効性の追加検証 |
| #93 | open | #10 context | flagellaのbody貫通検証 |

MVPでは #124 / #69 / #15 は Phase 3 baseline をblockしない。必要なら別系統で継続する。

### Early / Historical Phase 2 Anchors

| Issue | State | Role |
| ---: | --- | --- |
| #5 | closed | 物理モデル実装方法の調査 |
| #7 | closed | 菌体・べん毛要件調査 |
| #16 | closed | 論文モデルの力学改善 |
| #19 | closed | べん毛数custom |
| #20 | closed | run/tumble状態の明示 |
| #21 | closed | 2D投影の菌体中心化 |
| #23 | closed | 旧べん毛束化条件検討 |
| #25 | closed | 菌体-第1べん毛ビーズ制約 |
| #27 / #29 / #33 | closed | projection削除・整理 |
| #36 | closed | 3D plot上のspring / hook表記 |
| #37 | closed | トルクあり条件での安定化 |
| #43 | closed | 滑らかな回転安定性 |
| #51 | closed | べん毛根元ねじれ由来torque伝搬 |
| #53 / #54 | closed | 単一べん毛 torque伝搬モデル評価・まとめ |
| #58 | closed | 後方束化による形状維持・遊泳検証 |

これらは通常作業では直接触らず，Phase 2の歴史的根拠が必要なときだけ `docs/phase2/phase2_tasks.md` と `docs/codex-runs/*/review_result.json` から読む。

### Tooling / Workflow

| Issue | State | Role |
| ---: | --- | --- |
| #35 | open | Codex開発環境整理 |
| #44 | closed | Codex開発・レビュー分離workflow |
| #50 | closed | SKILL.md追加 |
| #64 | closed | CodexによるPRレビュー |
| #96 | closed | `scripts/01_simulate_swimming` 整理 |
| #101 | closed | pre-commit軽量化 |
| #106 | closed | multi-run / sweep / heatmap / replay共通導線 |
| #109 | open | plot / replay CLI に run.log / manifest追加 |
| #110 | open | replay CLI の出力内容選択 |

## User-Run Required Queue

現時点でユーザー実行が必須の項目はない。

今後ユーザー実行が必要になった場合は，この形式で報告する。

```bash
uv run python ...
```

確認先:

```text
outputs/...
```

判断ポイント:

```text
body deformation が見えるか
helical shape が保たれているか
2D pseudo-microscopy として自然か
同一run由来clipを独立sampleとして扱っていないか
```

## Next Three Actions

1. #127 で `source_video` / `track` / `clip` / `frame` のschemaとgroup keyを固定する。
2. #129 で `0.25 s`, `0.5 s`, `1.0 s` clip の評価計画と独立run数の数え方を決める。
3. #126 の初期2D投影特徴量解析を受けて，frame由来特徴量で同じ傾向が抽出できるか評価する。
