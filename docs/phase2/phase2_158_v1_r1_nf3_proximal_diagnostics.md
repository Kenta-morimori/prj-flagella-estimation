# Phase 2 Issue #158: dataset v1 r1 n=3 proximal failure diagnosis

## Scope

この文書は dataset v1 r1 の 3.0 s RUN 固定 raw output を再解析し，`n_flagella=3` で発生した非定常回転と proximal flag bond failure の原因候補を整理する。

今回の対象は原因解明と解決策の整理までであり，dataset v2 の生成，5 s x 27 run quality gate，Phase 4 ML MVP 実装は行わない。

入力:

- raw campaign: `outputs/phase2_multi_run/flagella_count_duration_3s_r1`
- campaign config: `conf/phase2_multi_run/flagella_count_duration_3s_r1.yaml`
- dataset_version: `v1`
- dataset_revision: `r1`
- source duration: `3.0 s`
- conditions: `n_flagella=1,2,3` x `attach_seed=0,1,2` x `phase_seed=0,1,2`

## Diagnostic CLI

既存の `step_summary.csv` / `state_archive.npz` を読み，長時間simulationは再実行しない。

```bash
uv run python scripts/03_dataset_building/analyze_dataset.py --analysis-kind diagnose-158 \
  --output-dir outputs/2026-07-31/phase2_158_probe \
  --overwrite
```

主な出力:

- `run_diagnostic_summary.csv`: 27 run の PASS/FAIL，first-fail，flag ID，local bond，seed別指標
- `failure_event_table.csv`: 4 fail条件の first-fail 前後同期表
- `failure_lead_lag_summary.csv`: first-fail 前 0.25 / 0.10 / 0.05 s の bond growth と contact proxy 先行性
- `seed_failure_table.csv`: `n_flagella=3` の attach / phase seed 別 failure 表
- `attach_geometry_table.csv`: attach body bead，attach pair distance，attach centroid offset
- `plots/<condition_id>_first_fail_sync.png`: angular velocity / speed / proximal bond / contact proxy の同期plot
- `hypothesis_assessment.json`
- `hypothesis_assessment.md`
- `manifest.json`
- `run.log`

replay は既存 `state_archive.npz` から生成する。長時間simulationは再実行しない。

```bash
uv run python scripts/03_dataset_building/replay_dataset.py --run-dir \
  --config conf/phase2_multi_run/flagella_count_duration_3s_r1.yaml \
  --input-dir outputs/phase2_multi_run/flagella_count_duration_3s_r1 \
  --output-dir outputs/2026-07-31/phase2_158_probe/replay \
  --mode both \
  --fps-out-3d 5 \
  --max-panels-per-grid 9 \
  --overwrite
```

replay output:

- `outputs/2026-07-31/phase2_158_probe/replay/grid_swim3d_page01.mp4`
- `outputs/2026-07-31/phase2_158_probe/replay/grid_swim3d_page02.mp4`
- `outputs/2026-07-31/phase2_158_probe/replay/grid_swim3d_page03.mp4`
- `outputs/2026-07-31/phase2_158_probe/replay/shape_stability_metrics.png`

## Result Summary

probe output: `outputs/2026-07-31/phase2_158_probe`

strict-pass:

| n_flagella | strict-pass |
|---:|---:|
| 1 | 9/9 |
| 2 | 9/9 |
| 3 | 5/9 |

`n_flagella=3` の fail は4条件で，いずれも `flag` category だった。

| condition_id | first_fail_t_s | first_fail_flag_id | first_fail_local_bond | max_flag_bond_rel_err | max local bond |
|---|---:|---:|---|---:|---|
| `as001__ps001__nf03` | 1.7803 | 1 | `2-3` | 1.2170 | `2-3` |
| `as002__ps000__nf03` | 2.6816 | 2 | `1-2` | 1.0265 | `1-2` |
| `as002__ps001__nf03` | 1.1870 | 1 | `1-2` | 1.2314 | `1-2` |
| `as002__ps002__nf03` | 2.4446 | 0 | `1-2` | 1.7206 | `1-2` |

`n_flagella=3` の attach geometry:

| attach_seed | attach body indices | min pair distance [um] | max pair distance [um] | centroid offset [um] | result pattern |
|---:|---|---:|---:|---:|---|
| 0 | `6|7|8` | 0.8660 | 0.8660 | ~0 | 3/3 pass |
| 1 | `6|7|3` | 0.5000 | 1.0000 | 0.3333 | 2/3 pass, 1/3 fail |
| 2 | `6|7|4` | 0.5000 | 1.0000 | 0.3333 | 0/3 pass, 3/3 fail |

first-fail 前 0.05 s の lead/lag summary:

| condition_id | bond delta | slope / s | flag-flag repulsion max | flag-body min gap | interpretation |
|---|---:|---:|---:|---:|---|
| `as001__ps001__nf03` | 0.1547 | 3.6819 | `3.69e-14` | 0.0284 | bond growth with contact correlation |
| `as002__ps000__nf03` | 0.4518 | 9.4610 | 0.0 | 0.0259 | bond growth without contact precursor |
| `as002__ps001__nf03` | 0.2296 | 3.6789 | 0.0 | -0.00017 | bond growth with contact correlation |
| `as002__ps002__nf03` | 0.5457 | 10.6414 | 0.0 | 0.0077 | bond growth without contact precursor |

## Proximal Relief Diagnostic

ユーザー実行の追加比較として，4 fail条件だけを `stiffness_scales.proximal_flag_spring=2.0` で3.0 s再実行した。これは dataset v2 生成ではなく，failure位置である flag chain local bond `1-2` / `2-3` だけを補強した原因分離probeである。

実行config:

- `conf/phase2_multi_run/flagella_count_v1_r1_nf3_proximal_relief_fail4.yaml`
- output: `outputs/phase2_multi_run/flagella_count_v1_r1_nf3_proximal_relief_fail4_pfs2_full`
- base condition: dataset v1 r1 3.0 s と同一。ただし `stiffness_scales.proximal_flag_spring=2.0`。

結果:

| condition_id | final_t_s | final pass | first_fail_t_s | max_flag_bond_rel_err | max local bond | body pass |
|---|---:|---|---:|---:|---|---|
| `as001__ps001__nf03__pfs2` | 2.9999 | True | none | 0.6067 | `4-5` | True |
| `as002__ps000__nf03__pfs2` | 2.9999 | True | none | 0.4153 | `1-2` | True |
| `as002__ps001__nf03__pfs2` | 2.9999 | True | none | 0.8695 | `3-4` | True |
| `as002__ps002__nf03__pfs2` | 2.9999 | True | none | 0.5114 | `3-4` | True |

baselineでは4条件すべてが `flag` category でfailし，first-fail local bondは `1-2` または `2-3` だった。一方，`proximal_flag_spring=2.0` では4条件すべてが3.0 sまで `final_shape_pass_nonbody=True`，`first_fail_t_s` なし，`body_shape_pass=True` になった。

したがって，proximal local flag bond `1-2` / `2-3` の局所剛性・load path が failure の最有力修正対象である。contact は一部条件で副次的な負荷変動に関与しうるが，今回の `pfs2` 結果は「接触回避」よりも「basal load transfer から proximal bond への局所変形集中」が主因であるという解釈を強める。

ただしユーザー定性評価では，`pfs2` でも菌体後方での自然な束化は確認できなかった。したがって `pfs2` は proximal flag bond failure の原因分離と形状破綻抑制には有効だが，dataset v2 の採用条件としては未確定である。dataset v2 では，他要件と合わせて滑らかな回転，泳ぎらしい後方束化，過剛性の副作用を別途評価する。

## Evidence

### Attach Seed

`attach_seed=0` は `n=3` の3 phaseすべて strict-pass した。`attach_seed=2` は3 phaseすべて fail し，`attach_seed=1` は `phase_seed=1` のみ fail した。

`attach_seed=0` は中心層の対称三角形配置で，attach pair distance は3組とも `0.866 um`，centroid offset はほぼ0だった。`attach_seed=1/2` は最小pair distance `0.5 um` の近接attach pairと centroid offset `0.333 um` を持つ非対称配置である。

したがって attach geometry は proximal bond 負荷集中の強い支配因子である。特に `attach_seed=2` は phase に依らず failure を誘発する。ただし `attach_seed=1` は phase依存でpass/failが分かれるため，非対称attach geometry単独ではなく，phase seed による初期helix位相と組み合わさって proximal bond への負荷集中が決まる。

### Phase Seed

同じ `attach_seed=1` でも `phase_seed=0,2` は pass，`phase_seed=1` は fail した。`attach_seed=2` の fail時刻と対象flagも phase により変わる。

したがって phase seed は failure の有無，時刻，対象flagを変える副因子である。

### Proximal Bond

4 failすべてで first-fail local bond は `1-2` または `2-3` だった。pass条件でも最大値は多くが proximal bond に出るが，fail条件では strict limit 付近を超える。

これは distal helix全体の崩壊ではなく，basal/root近傍の load transfer が proximal flag bond に集中していることを示す。

### Contact Proxy

`flag_flag_repulsion_force_max_N` は first-fail 前 0.05 s で，4 fail中3条件が `0.0` だった。`as001__ps001__nf03` だけは `3.69e-14 N` 程度の flag-flag repulsion がある。

flag-body は bead center distance から `2a` を引いた proxy gap で確認した。小さな負値は pass条件にも見えるため，現時点では flag-body contact / penetration を failure の主因とは判定しない。ただし bodyを連続capsuleではなくbead proxyで見ているため，完全な貫通判定ではない。

0.05 s lead/lag では，4 fail中2条件が `bond_growth_without_contact_precursor`，2条件が `bond_growth_with_contact_correlation` だった。したがって contact は一部条件の副次的な負荷変動には寄与しうるが，4 fail全体を説明する主因ではない。

### Body Roll And Axis Spin

`axis_center_to_body_roll_ratio_mean` は `n=3` で大きくseed依存する。fail条件だけに一意ではないが，`as001__ps001__nf03` は `57.6` と特に高い。非定常回転は proximal bond failure と同時に読むべきで，body roll単独では failure を説明しない。

## Hypothesis Classification

| hypothesis | classification | reason |
|---|---|---|
| attach geometry が proximal bond へ局所負荷を集中させる | support | failが `attach_seed=1,2` に偏り，`attach_seed=2` は3/3 fail |
| phase seed が failure時刻・対象flag・回転変動を変える | support | 同一attach内でpass/failと時刻が変わる |
| flag-flag / flag-body contact が主因 | weak as primary cause | 0.05 s lead/lagで2/4 failはcontact precursorなし，flag-flag repulsionは3/4 failで0 |
| attach-frame拘束，トルク伝達，body roll，axis spin の不均衡 | support | seed依存の回転比とproximal負荷が同期しており，basal load transferの問題と整合 |
| 物理的変動とartifactが混在 | support | attach/phase依存は物理的ばらつき，proximal過伸長は現行拘束・数値モデルの限界 |

## Most Likely Cause

最有力原因は，`seeded_surface` の非対称attach geometry，初期helix phase，`root_torque_segment_couples` のbasal load transferが組み合わさり，`n=3` の一部配置で proximal local flag bond `1-2` / `2-3` に局所負荷を集中させることである。`proximal_flag_spring=2.0` の4 fail条件再実行で3.0 s first-failなしになったため，この proximal load path が第一修正対象である。

より具体的には，対称attach (`attach_seed=0`) では3 phaseすべて pass する一方，centroid offset `0.333 um` と近接attach pair `0.5 um` を持つ非対称attachでは proximal bond margin が小さくなり，`attach_seed=2` では全phaseで strict limit を超える。phase seed は対象flagとfailure時刻を変える。

flag-flag contact は一部条件で寄与しうるが，4 fail全体を説明する主因ではない。flag-body contact proxyも pass/failを分ける強い証拠にはなっていない。

## Solution Direction

dataset v2 生成前に検討すべき修正対象:

1. proximal flag bond 専用の load relief または stiffness profile を追加する。
   - 対象は hook `first-second` ではなく flag chain local `1-2` / `2-3`。
   - 既存の `motor.local_first_second_spring_scale` は hook/first-second近傍の別補強であり，今回の主failureを直接狙っていない。
2. attach geometry quality gate を追加する。
   - `attach_seed=2` のように phase非依存でfailする配置を，dataset生成前のgeometry QCで検出する。
   - ただし単純なseed除外はdataset varietyを狭めるため，最終策ではなく診断用 fallback とする。
3. torque distribution / basal freedom の変更は第二候補にする。
   - body roll / axis spin imbalance は見えるが，first-fail local bondが一貫してproximalなので，まず局所load pathを直す方が狭い修正になる。

今回の `pfs2` は「壊れない」ことを示す診断候補であり，「菌体後方で束化する」ことまでは満たしていない。dataset v2 で同様の形状破綻が見られた場合は，今回の知見を proximal load path の第一切り分け材料として使う。ただし新規Issueは作らず，dataset v2 の物理条件検討の中で扱う。

物理解釈を変えるため，1 の実装修正前にはユーザー確認が必要である。

## Next Minimal Comparisons

長時間27 runを再実行する前に，4 fail条件だけで以下を比較する。

1. baseline replay-only diagnostics:
   - existing outputを本CLIで再解析する。
2. proximal local flag spring candidate:
   - `local 1-2 / 2-3` の flag springだけを段階的に `1.25`, `1.5`, `2.0` へ上げる診断extensionを追加する。
   - `2.0` はユーザー実行で4 fail条件すべて3.0 s passを確認済み。
   - 次に dataset v2 要件と合わせる前に，`1.25` / `1.5` の最小倍率確認，`n=1,2,3` 全27条件のsmoke quality gate，滑らかな回転と後方束化の定性評価を後続作業で扱う。
3. attach geometry reject-only reference:
   - モデル修正なしで `attach_seed=2` を除いた場合にdataset coverageがどう落ちるかを表にする。

上記のうち simulation を伴う比較は長時間run扱いになるため，Codexはユーザーの明示指示なしには実行しない。
