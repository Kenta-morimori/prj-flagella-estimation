# Phase 3.3 dataset mixing and versioning rules

この文書は Issue #128 の設計メモである。Phase 2 由来の条件変更を，観測 augmentation，独立 run，domain variation，dataset version 変更へ分類し，#127 schema の provenance と split group に接続する。

## Scope

ここで決める:

- 条件変更の4分類
- mixed training dataset に入れてよい条件と入れない条件
- `model_id`, `render_id`, `dataset_version`, `group_key` の責務
- torque / Brownian / RUN-TUMBLE / render variation / model変更の初期扱い
- MVP / v2 の固定 decision point

ここで決めない:

- 実顕微鏡動画を training label 付き dataset として採用するか
- Phase 2 物理モデルの変更

## Fixed Decisions

- MVP training baseline は Phase 2 dataset v1 の RUN固定 `n_flagella=1,2,3` を維持する。
- torque variation は MVP training baseline に混ぜず，diagnostic / robustness dataset として分離する。
- Brownian は当面含めない。現時点では `n_flagella=0` や thermal fluctuation の検証に効く条件であり，MVP の `n_flagella=1,2,3` 範囲外とする。
- RUN-TUMBLE は含めたいが，v2以降の別 dataset とする。
- render variation は軽い観測 augmentation のみ training に含める。広い domain adaptation は robustness / validation 側に分離する。
- `n_flagella=4` は v1/MVP では diagnostic-only のまま，v2 で再検討する。

## Diagnostic-only

`diagnostic-only` は，学習用 train / validation / test split には入れず，問題の切り分け，将来改善の根拠，robustness 確認だけに使う条件を指す。metadata や summary は保存してよいが，Phase 4 の正式 training candidate として扱ってはいけない。

代表例は v1 の `n_flagella=4` である。9 sample 中6 sample は strict pass したが，3 sample に `flag` failure と非等速回転が残ったため，v1/MVP training label には入れない。v2 では #124 の物理改善や RUN-TUMBLE profile の整理後に再検討する。

## Classification

| class | examples | label-preserving? | split grouping | versioning |
| --- | --- | --- | --- | --- |
| 観測 augmentation | crop jitter, rotation, translation, brightness, contrast, noise, blur, codec, render profile | 条件範囲内なら yes | base `group_key` を維持 | 原則 `dataset_version` は上げない。`render_id` / augmentation metadata に記録 |
| 独立 run | attach seed, phase seed, initial body pose, independent simulation seed | yes if same model and behavior regime | run ごとに別 `group_key` | 原則 `dataset_version` は上げない |
| domain variation | torque, viscosity, temperature proxy, Brownian, RUN / TUMBLE, frame-rate policy that changes motion evidence | conditional | variation 元 run 単位で group。base run 派生なら同 group | 採用範囲を registry に記録。training baseline に混ぜる場合は dataset freeze checklist が必要 |
| dataset version 変更 | equations, force distribution, constraint, stiffness profile, geometry, body/basal freedom model, label range, training candidate range | not automatically comparable | version をまたいだ leakage も audit | `dataset_version` と通常 `model_id` を更新 |

augmentation は観測過程の変換であり，細胞運動そのものを変える条件ではない。独立 run は同じ分布からの追加観測である。domain variation は現実的範囲と label-preserving 性の根拠がない限り，通常 training dataset へ混ぜない。

## Schema Connections

| field | role | rule |
| --- | --- | --- |
| `provenance.model_id` | 物理モデル・主要 stiffness・force distribution など dataset 解釈に影響する条件 | model変更時に更新する。同じ training dataset 内で複数 model を混ぜる場合は明示的な承認が必要 |
| `provenance.dataset_version` | dataset 論理版 | Phase 2 registry と一致させる。model変更，label range変更，baseline採用範囲変更で更新 |
| `provenance.render_id` | 同一 raw run 由来の描画条件 | render variation ごとに変えてよいが，split group には使わない |
| `track.group_key` | split leakage 防止の正本 | 同一 raw run / source track / render variation / overlap window は同じ group に置く |

推奨 `group_key` は #129 と同じく pseudo では `phase2:<dataset_version>:<run_id>`，実動画では `real:<source_video_id>:<track_id>` とする。

## Mixing Rules

同じ baseline training dataset へ混ぜてよい初期条件:

- 同一 `dataset_version` / `model_id` の attach seed と phase seed 違い。
- 同一 raw run からの label-preserving な軽い観測 augmentation。ただし同じ `group_key` に束ねる。
- #8 / #17 の実動画条件に根拠を持つ小範囲の render variation。blur / noise / brightness / contrast など，2D特徴量を不自然に消さない範囲に限る。
- 欠損や低品質 clip を `qc.status` / `exclusion_reason` で除外できる範囲の処理差分。

混ぜない条件:

- torque 数値だけを変えた run を，現実的範囲とラベル保持の根拠なしに augmentation として混ぜること。
- Brownian あり / なしを同じ baseline として無条件に混ぜること。
- RUN 固定と RUN-TUMBLE を同じ class distribution として無条件に混ぜること。
- force distribution, constraint, stiffness profile, body geometry, basal freedom model の変更を同じ `dataset_version` として混ぜること。
- `n_flagella=4` diagnostic-only や `n_flagella>=5` を MVP training candidate へ黙って戻すこと。
- 同一 raw run 由来の render variation や overlap window を train / validation / test にまたがせること。

条件付きで混ぜられる候補:

- torque / viscosity の小範囲 variation: MVP では training に混ぜない。実測または文献上の想定範囲，2D特徴量の class separability，failure率，visual review を満たした後，diagnostic / robustness dataset として評価する。
- Brownian: 当面は training / v2 初期ともに含めない。#15 の実装・validation 後，thermal noise が対象タスクに必要になった場合だけ別 dataset として扱う。
- RUN-TUMBLE: #69 の実装後，v2以降の別 dataset として扱う。MVP baseline へは混ぜない。
- render variation: #17 の render profile が実動画条件に基づき，2D特徴量を不自然に消さないことを確認した小範囲だけ MVP training に含める。

## RUN-TUMBLE v2 Profile

論文例 config は `run_tumble.run_tau=1200`, `tumble_tau=800`, `semicoiled_tau=400`, `curly1_tau=400` である。現行 Phase 2 実装では `tau_s=1.0 s` のため，この値をそのまま読むと次の時間が必要になる。

| event | required duration |
| --- | ---: |
| first TUMBLE start | `>1200 s` |
| semicoiled phase end | `1600 s` |
| curly1 phase end / one full cycle | `2000 s` |

これは Phase 3 / 4 dataset 生成には長すぎる。v2 では論文の状態遷移構造，つまり RUN -> motor reversal -> semicoiled -> curly1 -> RUN を使いつつ，数秒以内で RUN / TUMBLE が観測できる短縮 profile を別定義する。短縮 profile は論文例そのものではなく dataset-generation profile として扱い，`dataset_version=v2` と `model_id` / `run_tumble_profile_id` 相当の provenance に記録する。

## Dataset Freeze Checklist

Phase 4 training dataset を凍結する前に確認する:

- `dataset_version`, `dataset_revision`, `dataset_id`, `model_id` が registry と一致している。
- `n_flagella` の training candidate 範囲が明記されている。
- 全 clip に `track.group_key` があり，grouped split leakage test が通る。
- `render_id` 違いと overlap window が split をまたいでいない。
- torque variation / RUN-TUMBLE / Brownian などの domain variation が MVP baseline へ混入していない。
- v2 以降で domain variation を含める場合，根拠 doc とユーザー承認がある。
- `qc.status=fail` の扱いと除外理由が固定されている。
- 実動画を含む場合，label source と利用目的が `unavailable`, `manual`, `phase2_gt` のどれかとして明確である。

## Remaining Decision Points

現時点で未決定として残す項目:

1. v2 の短縮 RUN-TUMBLE profile の具体的な `run_tau`, `tumble_tau`, `semicoiled_tau`, `curly1_tau`。
2. v2 で `n_flagella=4` を再検討するとき，#124 の物理改善完了を必須にするか。
3. render variation の数値範囲。#8 / #17 の実動画条件を根拠に決める。

MVP baseline では Phase 2 dataset v1 の RUN固定 `n_flagella=1,2,3`，同一 model の独立 run と小範囲の render augmentation だけを混ぜる。torque / Brownian / RUN-TUMBLE / model変更は別 dataset または diagnostic として扱う。
