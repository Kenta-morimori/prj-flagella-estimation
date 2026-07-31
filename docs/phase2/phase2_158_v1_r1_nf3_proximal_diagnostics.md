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
uv run python scripts/02_phase2_analysis/diagnose_v1_r1_nf3_failures.py \
  --output-dir outputs/2026-07-31/phase2_158_probe \
  --overwrite
```

主な出力:

- `run_diagnostic_summary.csv`: 27 run の PASS/FAIL，first-fail，flag ID，local bond，seed別指標
- `failure_event_table.csv`: 4 fail条件の first-fail 前後同期表
- `plots/<condition_id>_first_fail_sync.png`: angular velocity / speed / proximal bond / contact proxy の同期plot
- `hypothesis_assessment.json`
- `hypothesis_assessment.md`
- `manifest.json`
- `run.log`

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

## Evidence

### Attach Seed

`attach_seed=0` は `n=3` の3 phaseすべて strict-pass した。`attach_seed=2` は3 phaseすべて fail し，`attach_seed=1` は `phase_seed=1` のみ fail した。

したがって attach geometry は proximal bond 負荷集中の強い支配因子である。特に `attach_seed=2` は phase に依らず failure を誘発する。

### Phase Seed

同じ `attach_seed=1` でも `phase_seed=0,2` は pass，`phase_seed=1` は fail した。`attach_seed=2` の fail時刻と対象flagも phase により変わる。

したがって phase seed は failure の有無，時刻，対象flagを変える副因子である。

### Proximal Bond

4 failすべてで first-fail local bond は `1-2` または `2-3` だった。pass条件でも最大値は多くが proximal bond に出るが，fail条件では strict limit 付近を超える。

これは distal helix全体の崩壊ではなく，basal/root近傍の load transfer が proximal flag bond に集中していることを示す。

### Contact Proxy

`flag_flag_repulsion_force_max_N` は first-fail 前 0.25 s で，4 fail中3条件が `0.0` だった。`as001__ps001__nf03` だけは `3.73e-14 N` 程度の flag-flag repulsion がある。

flag-body は bead center distance から `2a` を引いた proxy gap で確認した。小さな負値は pass条件にも見えるため，現時点では flag-body contact / penetration を failure の主因とは判定しない。ただし bodyを連続capsuleではなくbead proxyで見ているため，完全な貫通判定ではない。

### Body Roll And Axis Spin

`axis_center_to_body_roll_ratio_mean` は `n=3` で大きくseed依存する。fail条件だけに一意ではないが，`as001__ps001__nf03` は `57.6` と特に高い。非定常回転は proximal bond failure と同時に読むべきで，body roll単独では failure を説明しない。

## Hypothesis Classification

| hypothesis | classification | reason |
|---|---|---|
| attach geometry が proximal bond へ局所負荷を集中させる | support | failが `attach_seed=1,2` に偏り，`attach_seed=2` は3/3 fail |
| phase seed が failure時刻・対象flag・回転変動を変える | support | 同一attach内でpass/failと時刻が変わる |
| flag-flag / flag-body contact が主因 | inconclusive, weak as primary cause | flag-flag repulsionは3/4 failでfirst-fail前0，flag-body proxyはpassにも出る |
| attach-frame拘束，トルク伝達，body roll，axis spin の不均衡 | support | seed依存の回転比とproximal負荷が同期しており，basal load transferの問題と整合 |
| 物理的変動とartifactが混在 | support | attach/phase依存は物理的ばらつき，proximal過伸長は現行拘束・数値モデルの限界 |

## Most Likely Cause

最有力原因は，`seeded_surface` attach geometry と `root_torque_segment_couples` のbasal load transferが組み合わさり，`n=3` の一部配置で proximal local flag bond `1-2` / `2-3` に局所負荷を集中させることである。

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

物理解釈を変えるため，1 の実装修正前にはユーザー確認が必要である。

## Next Minimal Comparisons

長時間27 runを再実行する前に，4 fail条件だけで以下を比較する。

1. baseline replay-only diagnostics:
   - existing outputを本CLIで再解析する。
2. proximal local flag spring candidate:
   - `local 1-2 / 2-3` の flag springだけを段階的に `1.25`, `1.5`, `2.0` へ上げる診断extensionを追加する。
   - 4 fail条件，3.0 s，strict QC，body roll / axis spin / speed特徴保持を確認する。
3. attach geometry reject-only reference:
   - モデル修正なしで `attach_seed=2` を除いた場合にdataset coverageがどう落ちるかを表にする。

上記のうち simulation を伴う比較は長時間run扱いになるため，Codexはユーザーの明示指示なしには実行しない。
