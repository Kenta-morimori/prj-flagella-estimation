# ADR 0009: Phase 2 2010 potential formulation と既存互換

- status: accepted
- date: 2026-08-02
- scope: Phase 2 / Issue #163

## Context

Watari & Larson (2010) の FENE-Fraenkel spring は、bond平衡長に対する相対変形と分母一乗を使う。現行springは `s*b` の絶対変形限界と分母二乗を使っており、式と次元が一致しない。一方、現行springを直接置換すると過去runを再現できない。

bending、torsion、hook potentialは論文表示の負号と現行実装の正号が異なる。ただし現行実装は平衡角への復元力として動作し、論文式の負号だけを通常の `F=-grad U` に適用すると平衡点から離す向きになる。

spring-spring repulsionは論文がexponential potentialとTable 1の値を与えるが、対象segmentと端点への力分配を完全には定めていない。

式・次元・実装対応の正本は `docs/phase2/phase2_163_2010_potential_correspondence.md` とする。

## Decision

- `potentials.spring.formulation` に `legacy` と `fene_fraenkel` を追加する。
- selector欠落時は `legacy` とし、過去runの数値挙動を維持する。
- `fene_fraenkel` は `q=(r-L)/L`、`|q|<s`、`H q/(1-q^2/s^2)` を使い、`H=10 T/b` をforceとして扱う。
- `legacy` は絶対限界 `s*b`、分母二乗、既存clipを変更しない。論文準拠ではなく互換用であることを明記する。
- bending、torsion、hookの符号は変更せず、平衡角への復元方向をtestで仕様化する。torsionの現行dihedral符号、`[-pi, pi)` wrap、有限差分も維持する。
- hookは `theta>90 deg` で無効、`theta<=90 deg` で90度へ復元する。
- repulsionは共有端点segmentを除外し、body-bodyを除外してbody-flagellum/flagellum-flagellumへ適用し、最短点の力を4端点へ線形分配する現行仮定を維持する。
- repulsionのcutoff不連続と完全交差時の方向退化は既知の数値近似として受容する。

## Default decision gate

比較結果がない段階では `legacy` をdefaultとする。motor-off `0.1 tau` と motor-on RUN `1 tau` の両方について、`fene_fraenkel` が完走し、非有限値・step欠損がなく、body/nonbody strict shape gateを全期間通過した場合だけ、論文一致を優先して2010年baselineの採用候補にする。いずれかを満たさない場合は `legacy` を維持する。

2026-08-02のuser-runでは、両formulationがmotor-off/onとも予定stepを完走し、body/nonbody strict shape gateを全期間通過した。corrected decision artifact `outputs/phase2_potential_comparison/2026-08-02/154149/default_decision.json` は `fene_fraenkel` を選択したため、`conf/sim_swim.yaml` のdefaultを `fene_fraenkel` とする。selector欠落時のparser fallbackは過去互換のため `legacy` のまま維持する。

初回解析 `outputs/phase2_potential_comparison/2026-08-02/153848` は、`step_summary.csv` がstep実行前時刻を記録する契約を考慮せず、正常な最終時刻 `duration-dt` を途中停止と誤判定したため採否根拠に用いない。判定器は連続step番号、観測dt、予定durationを照合するよう修正した。

## Consequences

- 既存configは変更なしで読み込め、過去runを再現できる。
- 論文準拠springを明示的に比較できる。
- `legacy` と `fene_fraenkel` の `H` は同じconfig値でも式中の次元的役割が異なるため、結果比較ではformulationを必ず記録する。
- potentialの表示符号より、実装された力の復元方向を契約として優先する。
- `sim_swim_2010.yaml` の作成とconfig再編はIssue #164、正式な `tau` duration指定はIssue #165で扱う。Issue #163では現行 `tau_s=1.0` 規約を比較profileに限定して用いる。

## Verification status

unit test、full pytest、force-extension artifact、比較profileのdry-run、motor-off/onの長時間比較、自動default判定まで完了した。動画目視は採否基準に含めない。
