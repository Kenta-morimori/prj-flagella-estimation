# Phase 2 sweep profiles

この directory には，Phase 2 の sweep / heatmap 用 profile YAML を置く。
実行コマンドの詳細は `scripts/README.md` を参照する。

## Canonical profile

`shape_stability_grid` 系の現役 user-facing 導線は次の 2 つを正本とする。

* `conf/phase2_sweeps/shape_stability_grid.yaml`
* `conf/phase2_sweeps/shape_stability_heatmap.yaml`

新しい実行例，README，issue / PR の説明は原則としてこの 2 つを使う。

`hook_overstretch.yaml` / `hook_overstretch_heatmap.yaml` は historical alias であり，
既存のメモ，過去 run，再現コマンドとの互換用に残す。削除前提では扱わない。

## `kind` と `mode`

`kind` は「どの処理プログラムを使うか」を表す。
`mode` は，その `kind` 実装の中で「どの条件パターンを使うか」を表す。

例:

```yaml
kind: shape_stability_grid
args:
  mode: attach-frame-grid
```

この profile は `shape_stability_grid` 実装を使い，
その中で `attach-frame-grid` の条件展開を選ぶ。

heatmap profile でも `mode` は重要である。
heatmap 側では，`summary.csv` のどの `mode` 行を読み，
どの列を x/y 軸に使うかを決める。

## `shape_stability_grid` の主な `mode`

| mode | 用途 |
| --- | --- |
| `preset` | 代表的な local mitigation 条件の比較 |
| `body-first-grid` | body attach から first bead までの距離・角度補強の比較 |
| `first-second-grid` | first bead から second bead までの距離補強の比較 |
| `attach-frame-grid` | attach frame の position / tangent scale の比較 |
| `torque-profile-grid` | torque distribution / profile 組み合わせの比較 |
| `basal-freedom-grid` | basal freedom 診断条件の比較 |
| `position-only-grid` | attach frame の position scale だけの比較 |

## Profile を増やすときの目安

既存 profile の既定値を一時的に変えるだけなら，新しい YAML を作らず
`KEY=VALUE` override を使う。

新しい YAML を作るのは，繰り返し使う代表条件，または
結果を後から参照する固定条件セットに限る。

heatmap 用 profile は sweep 条件ごとに増やさず，
同じ heatmap 実装で読める場合は `summary_csv=...` と `mode=...` を実行時指定する。
