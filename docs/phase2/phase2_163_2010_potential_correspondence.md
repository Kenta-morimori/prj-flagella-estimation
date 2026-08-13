# Phase 2 Issue #163: Watari & Larson (2010) potential 式照合

## 目的

Watari & Larson (2010) の bead-spring model と現行実装を、式、変数、SI次元、復元方向まで照合する。論文式をそのまま転記した箇所、現行実装、project-specific extension、数値安定化近似を区別する。

参照一次資料:

- Nobuhiko Watari and Ronald G. Larson, "The hydrodynamics of a run-and-tumble bacterium propelled by polymorphic helical flagella", *Biophysical Journal* 98, 12-17 (2010), DOI: [10.1016/j.bpj.2009.09.044](https://doi.org/10.1016/j.bpj.2009.09.044)
- 論文本文の Model 節、Eqs. (1)-(4)、Table 1

## 記号と次元

論文では長さを `b`、時間を `tau = eta b^3 / T`、エネルギーを motor torque の大きさ `T` でscaleする。SIでは `b, L, r, D, alpha_ss` は `[m]`、`T, k_b, k_t, A_ss` は `[J] = [N m]`、`H = 10 T/b` は `[N]`、`s` は無次元である。角度は無次元量として扱う。

`r_ij = r_i - r_j`、`r_ij = |r_ij|`、`L` はbondごとの平衡長、`q = r_ij/L - 1` は相対変形とする。現行configの `H_over_T_over_b`, `kb_over_T`, `kt_over_T`, `A_ss_over_T`, `a_ss_over_b`, `cutoff_over_b` は各論文scaleに対する無次元比である。

## Spring

論文 Eq. (1) の FENE-Fraenkel force は、向きの規約を除く大きさとして次である。

```text
q = r_ij / L - 1
|F_FF| = H q / (1 - q^2 / s^2)
1 - s < r_ij / L < 1 + s
```

Table 1 は `H = 10 T/b`, `s = 0.1`、flagellum の `L = 0.58 b` を与える。したがって `s` は絶対長ではなく、各bondの平衡長に対する相対変形限界である。伸長時はbondを縮め、圧縮時はbondを伸ばす向きを復元方向として固定する。論文のベクトル式は `F_FF(r_ij)` がどちらのbeadに作用する力かを本文中で明示していないため、式先頭の符号だけでは現行配列への加算符号を決めない。

Issue #163 では `potentials.spring.formulation` を追加する。

| formulation | force magnitude | limit | 位置づけ |
| --- | --- | --- | --- |
| `fene_fraenkel` | `H q / (1 - q^2/s^2)` | `|q| < s` | 論文 Eq. (1) 準拠。`H` は `[N]` |
| `legacy` | `H_impl x / (1 - x^2/(s b)^2)^2` | `|x| < s b` | 過去run再現用。`x=r-L`、絶対限界、分母二乗を維持 |

`legacy` の式は論文の FENE-Fraenkel force と次元・関数形が一致しない。既存結果との数値互換だけを目的として変更しない。両式とも特異点で直接発散させず、許容端の内側へclipする。これは数値安定化近似であり、論文由来ではない。

既存configでselectorが欠落した場合と `formulation=legacy` を明示した場合は同一挙動とする。`H <= 0`、`s <= 0`、未知formulationを拒否し、`fene_fraenkel` では圧縮側の許容長を正に保つため `s < 1` も要求する。

## Bending, torsion, hook

論文 Eqs. (2), (3) は次の負号付きpotentialを表示する。

```text
phi_b(theta) = -0.5 k_b (cos(theta) - cos(theta_0))^2
phi_t(phi)   = -0.5 k_t (phi - phi_0)^2
```

一方、現行実装は正の二乗energyを用い、`F = -grad U` により平衡角へ戻す。

```text
U_b(theta) = +0.5 k_b (cos(theta) - cos(theta_0))^2
U_t(phi)   = +0.5 k_t wrap(phi - phi_0)^2
```

論文の負号を通常のpotential energyとしてそのまま `F=-grad(phi)` へ入れると、平衡角は極小ではなく極大になる。論文側のpotential/force符号規約だけでは現行符号の反転を正当化できないため、Issue #163 では現行の復元方向を単体テストで固定し、符号は変更しない。

torsion は4 beadのdihedral angleを現行定義で計算し、`phi-phi_0` を `[-pi, pi)` にwrapした正の二乗energyの勾配を中心有限差分で求める。有限差分幅 `fd_eps_over_b` は実装上の数値近似であり、論文 Table 1 のparameterではない。

hook は論文 Eq. (4) に従い、`theta_hook > 90 deg` ではpotentialを無効、`theta_hook <= 90 deg` では90度へ戻す。論文は有効域を `-0.5 k_b cos^2(theta_hook)` と書くが、現行実装は復元方向を満たす `+0.5 k_b cos^2(theta_hook)` 相当を用いる。境界 `90 deg` は有効側だが力はゼロになる。

Table 1 の値は `k_b=20 T`, `k_t=10 T`、normal/semicoiled/curly1 の平衡角がそれぞれ bending `142/90/105 deg`、torsion `-60/65/120 deg` である。現行configの値はこれと一致する。

## Spring-spring repulsion

論文本文は次のpotentialだけを規定する。

```text
phi_ss(D) = A_ss exp(-D/alpha_ss)
|F_ss| = (A_ss/alpha_ss) exp(-D/alpha_ss)
```

`D` は2 spring segment間の最短距離である。Table 1 の論文由来値は `A_ss=1.0 T`, `alpha_ss=0.2 b`, cutoff `0.2 b` である。cutoffでpotential/forceをshiftしていないため、現行力はcutoffで不連続になる。

論文はsegment列挙とbeadへの分配を完全には規定していない。現行実装では次を実装上の仮定として固定する。

- 端点を共有する隣接segmentは除外する。
- body-body segment対は除外し、body-flagellum と flagellum-flagellum に適用する。
- segment上の最短点parameterを使い、反対向きの力を各segmentの2端点へ線形分配する。これにより合力はゼロになる。
- 最短点が完全に一致する交差では方向ベクトルが退化し、現行実装は有限な反発方向を生成できない。この点はhard non-crossing constraintではなく数値近似として扱う。

## 対応表

| 対象 | 論文 | 現行/Issue #163後 | 判定 |
| --- | --- | --- | --- |
| spring `H` | `10 T/b`、force `[N]` | `fene_fraenkel`で同じscale | 論文準拠 |
| spring `s` | 相対変形限界 `0.1` | `fene_fraenkel`は相対、`legacy`は絶対 `s b` | selectorで分離 |
| spring 分母 | 一乗 | `fene_fraenkel`一乗、`legacy`二乗 | selectorで分離 |
| bending | 負号付きcos差二乗 | 正号付きcos差二乗 | 復元方向を優先して維持 |
| torsion | 負号付き角度差二乗 | 正号、wrap、有限差分 | 復元方向を優先、数値近似あり |
| hook | `>90 deg`無効、`<=90 deg`負号付きcos二乗 | 同じ有効域、正号相当 | 復元方向を優先して維持 |
| repulsion式/値 | exponential、`0.2 b` range/cutoff | 一致 | 論文準拠 |
| repulsion対象/分配 | 詳細規定なし | 共有端点除外、body-body除外、線形端点分配 | 実装仮定 |

## Default採否と比較run

`fene_fraenkel` のdefault採否は、force-extension比較だけでは決めず、同一条件の motor-off `0.1 tau` と motor-on RUN `1 tau` で判断する。Issue #163 では `tau_s=1.0` の現行legacy時間規約により、比較profileをそれぞれ `0.1 s`, `1.0 s` とする。汎用的な `tau` 指定は Issue #165 に分離する。

自動判定は次で固定する。

1. `fene_fraenkel` の両runが例外、非有限値、予定step欠損なしで完走する。
2. 両runが既存のbody/nonbody strict shape gateを全期間通過する。
3. 1と2を満たせば、論文一致を優先して `fene_fraenkel` を2010年baselineの採用候補とする。それ以外は `legacy` を維持する。
4. 両formulationが失敗した場合も、過去互換を優先して `legacy` を維持し、失敗条件を診断結果として保存する。

campaignがcondition例外で `summary.csv` 作成前に停止した場合も、解析CLIは欠損summaryをrun失敗として扱い、理由付きのlegacy判定artifactを生成する。

長時間比較はリポジトリ方針に従いユーザーが実行する。実装側の短時間test、dry-run、force-extension artifact検証が完了した後に、実行コマンド、出力先、確認ファイル、判定結果の取込コマンドを一括提示する。自動基準で採否するため、この比較自体に動画目視は要求しない。

実行導線は次の2 profileと解析CLIに固定する。

```bash
uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/spring_formulation_motor_off.yaml

uv run python scripts/01_simulate_swimming/run_multi_run.py \
  config=conf/phase2_multi_run/spring_formulation_motor_on.yaml

uv run python scripts/02_phase2_analysis/analyze_spring_formulations.py \
  motor_off_summary=<motor-off-run-root>/summary.csv \
  motor_on_summary=<motor-on-run-root>/summary.csv
```

run rootはそれぞれ `outputs/phase2_multi_run/spring_formulation_motor_off/` と `outputs/phase2_multi_run/spring_formulation_motor_on/` のJST timestamp配下に作られる。解析出力はデフォルトで `outputs/phase2_potential_comparison/` のJST timestamp配下に作られ、`force_extension.csv`, `force_extension.png`, `default_decision.json`, `default_decision.md`, `run.log`, `manifest.json` を確認する。

2026-08-02の採否runは次である。

- motor-off: `outputs/phase2_multi_run/spring_formulation_motor_off/2026-08-02/151207`
- motor-on: `outputs/phase2_multi_run/spring_formulation_motor_on/2026-08-02/151505`
- corrected decision: `outputs/phase2_potential_comparison/2026-08-02/154149`

両formulationは両runで完走し、body/nonbody strict shape gateを全期間通過した。`fene_fraenkel` は `legacy` に対して、motor-off/onの `max_flag_bond_rel_err` をそれぞれ49.2%/84.2%、`body_spring_max_stretch_ratio` を84.4%/81.0%低減した。採否規則と論文一致に従い、`conf/sim_swim.yaml` のdefaultを `fene_fraenkel` とする。selectorを省略した既存configは引き続き `legacy` へfallbackする。Issue #164 が扱う `sim_swim_2010.yaml` の作成・config再編は本Issueに含めない。

## Issue境界

- Issue #163: potential式照合、spring selector、復元方向test、repulsion仮定、比較導線と採否基準。
- Issue #164: `sim_swim_2010.yaml` を含む2010年baseline configの再編・命名。
- Issue #165: `tau` を明示したduration指定と時間scaleの正式なCLI/config契約。
