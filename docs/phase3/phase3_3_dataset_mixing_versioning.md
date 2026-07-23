# Phase 3.3 dataset mixing and versioning rules

この文書は Issue #128 の設計メモである。Phase 2 由来の条件変更を，観測 augmentation，独立 run，domain variation，dataset version 変更へ分類し，#127 schema の provenance と split group に接続する。

## Scope

ここで決める:

- 条件変更の4分類
- mixed training dataset に入れてよい条件と入れない条件
- `model_id`, `render_id`, `dataset_version`, `group_key` の責務
- torque / Brownian / RUN-TUMBLE / render variation / model変更の初期扱い
- ユーザー判断待ちの decision point

ここで決めない:

- torque variation などの採用可否の最終判断
- 実顕微鏡動画を training label 付き dataset として採用するか
- Phase 2 物理モデルの変更

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
- 同一 raw run からの label-preserving な観測 augmentation。ただし同じ `group_key` に束ねる。
- #8 / #17 の実動画条件に根拠を持つ render variation。split は raw run 単位で固定する。
- 欠損や低品質 clip を `qc.status` / `exclusion_reason` で除外できる範囲の処理差分。

混ぜない条件:

- torque 数値だけを変えた run を，現実的範囲とラベル保持の根拠なしに augmentation として混ぜること。
- Brownian あり / なしを同じ baseline として無条件に混ぜること。
- RUN 固定と RUN-TUMBLE を同じ class distribution として無条件に混ぜること。
- force distribution, constraint, stiffness profile, body geometry, basal freedom model の変更を同じ `dataset_version` として混ぜること。
- `n_flagella=4` diagnostic-only や `n_flagella>=5` を MVP training candidate へ黙って戻すこと。
- 同一 raw run 由来の render variation や overlap window を train / validation / test にまたがせること。

条件付きで混ぜられる候補:

- torque / viscosity の小範囲 variation: 実測または文献上の想定範囲，2D特徴量の class separability，failure率，visual review を満たす場合だけ domain variation として採用する。
- Brownian: #15 の実装・validation 後，thermal noise が実動画観測条件に必要で，label-preserving と判断できる場合だけ採用する。
- RUN-TUMBLE: #69 の実装後，Phase 4 の対象タスクが RUN 状態だけか RUN/TUMBLE 混在かをユーザーが決めた後に別 dataset として扱う。
- render variation: #17 の render profile が実動画条件に基づき，2D特徴量を不自然に消さないことを確認した範囲だけ採用する。

## Dataset Freeze Checklist

Phase 4 training dataset を凍結する前に確認する:

- `dataset_version`, `dataset_revision`, `dataset_id`, `model_id` が registry と一致している。
- `n_flagella` の training candidate 範囲が明記されている。
- 全 clip に `track.group_key` があり，grouped split leakage test が通る。
- `render_id` 違いと overlap window が split をまたいでいない。
- domain variation を含める場合，根拠 doc とユーザー承認がある。
- `qc.status=fail` の扱いと除外理由が固定されている。
- 実動画を含む場合，label source と利用目的が `unavailable`, `manual`, `phase2_gt` のどれかとして明確である。

## User Decision Points

現時点でユーザー判断待ちにする項目:

1. torque variation を Phase 4 training baseline に混ぜるか，別 diagnostic dataset に留めるか。
2. Brownian / RUN-TUMBLE を MVP baseline に含めるか，Phase 2 physics extension の成果後まで分離するか。
3. render variation を実動画 domain adaptation として training に含めるか，validation-only / robustness check に留めるか。
4. `n_flagella=4` を diagnostic-only のまま維持するか，#124 の物理改善後に別 version で再検討するか。

初期推奨は，MVP baseline では Phase 2 dataset v1 の RUN固定 `n_flagella=1,2,3`，同一 model の独立 run と根拠付き render augmentation だけを混ぜ，torque / Brownian / RUN-TUMBLE / model変更は別 dataset または diagnostic として扱うことである。
