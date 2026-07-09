# Phase 2 sweep profiles

この directory には，Phase 2 の sweep / heatmap 用 profile YAML を置く。
実行コマンドの詳細は `scripts/README.md` を参照する。

## Profile の基本形

```yaml
kind: shape_stability_grid
args:
  config: conf/sim_swim.yaml
  mode: preset
  output_dir: outputs/phase2_sweeps/shape_stability_grid
```

`kind` は，この profile をどの種類の処理として実行するかを表す。
たとえば `shape_stability_grid` は shape stability 用の sweep / heatmap 実装を使い，
`motor_scale` は motor scale 用の sweep 実装を使う。

`kind` が変わると，呼ばれる Python 実装が変わる。そのため，使える `args`，
出力される `summary.csv` の列，heatmap の作り方，評価対象も変わる。

`args` は，選ばれた `kind` の実装へ渡す引数である。
profile 内では YAML として書き，実行時には `--mode preset` のような CLI 引数へ変換される。
コマンドラインの `KEY=VALUE` は profile の既定値を上書きするために使う。

## `kind` と `mode`

`kind` は「どの処理プログラムを使うか」である。
`mode` は，その処理プログラムの中で「どの条件パターンを使うか」である。

`mode` は全 `kind` 共通の変数ではない。
同じ `mode` 名でも，別の `kind` では意味を持たないか，そもそも受け付けられない。

例:

```yaml
kind: shape_stability_grid
args:
  mode: attach-frame-grid
```

この profile は `shape_stability_grid` 実装を使い，その中で
`attach-frame-grid` という条件パターンを選ぶ。

heatmap profile でも `mode` は重要である。
heatmap 側では，`summary.csv` のどの `mode` 行を読むか，どの列を x/y 軸に使うかを決める。
sweep 側に存在する `mode` でも，heatmap 実装側が未対応なら heatmap では使えない。

## 主な `kind`

| kind | 何をするか |
| --- | --- |
| `motor_scale` | motor-local scale を変えた sweep を実行する |
| `single_flagellum_torque` | single flagellum の torque 条件を評価する |
| `bundling_alignment` | 複数べん毛の helix axis alignment を評価する |
| `shape_stability_grid` | hook / proximal flagellum を含む shape stability 条件を評価する |
| `motor_scale_collapse` | motor scale sweep の collapse heatmap を作る |
| `dt_star_torque` | `dt_star` x torque heatmap を作る |
| `local_scale_mode` | local scale mode x torque heatmap を作る |

## `shape_stability_grid` の主な `mode`

| mode | 何が変わるか |
| --- | --- |
| `preset` | 代表的な local mitigation 条件をまとめて比較する |
| `body-first-grid` | body attach から first bead までの距離・角度補強を振る |
| `first-second-grid` | first bead から second bead までの距離補強を振る |
| `attach-frame-grid` | attach frame の position / tangent scale を振る |
| `torque-profile-grid` | torque distribution / profile の組み合わせを比較する |
| `basal-freedom-grid` | basal freedom 診断用に position / tangent / torque 条件を比較する |
| `position-only-grid` | attach frame の position scale だけを振る |

## Profile を増やすときの目安

既存 profile の値だけを一時的に変える場合は，新しい YAML を作らず
`KEY=VALUE` override を使う。

新しい YAML を作るのは，繰り返し使う代表条件，または結果を後から参照する必要がある
固定条件セットに限る。

heatmap 用 profile は，原則として sweep 条件ごとに増やさない。
同じ heatmap 実装で読める場合は，`summary_csv=...` と `mode=...` を実行時に指定する。
