# Phase 2 実行コマンド リファレンス

本ドキュメントは、Phase 2 における主要な実行シーケンスの統一的なコマンド形式を提供する。

---

## 前提

- 全コマンドは workspace root で実行
- Python 環境: `uv` を使用（`.venv` 自動アクティベーション）
- 出力先: `outputs/` 配下（タイムスタンプ自動生成、または明示指定）
- 詳細なパラメータ説明は[prompts/phase2/0037.md](../../prompts/phase2/0037.md)を参照

---

## セクション 1: 単一実行（01_simulate_swimming.py）

### 基本構文

```bash
uv run python -m scripts.01_simulate_swimming \
  [--duration-s DURATION] \
  [--render-flagella] \
  [--no-render-flagella] \
  [KEY=VALUE]...
```

### パラメータ説明

| オプション | 型 | デフォルト | 説明 |
|-----------|-----|----------|------|
| `--duration-s` | float | 0.1 | シミュレーション時間（秒） |
| `--render-flagella` | flag | False | べん毛 3D レンダリング有効化 |
| `--no-render-flagella` | flag | - | べん毛 3D レンダリング無効化 |
| `KEY=VALUE` | 任意 | - | Hydra オーバーライド（`conf/sim_swim.yaml` の任意 key を上書き） |

### よく使う単一実行例

#### PhaseA 検証: Torque ゼロ，minimal_basal_stub (0.01s)

```bash
uv run python -m scripts.01_simulate_swimming \
  --duration-s 0.01 \
  --no-render-flagella \
  flagella.n_flagella=1 \
  flagella.stub_mode=minimal_basal_stub \
  motor.torque_Nm=0.0 \
  output.base_dir=outputs/verify_phaseA_torque0_minimal_1beads
```

出力先: `outputs/{タイムスタンプ}/` また `outputs/verify_phaseA_torque0_minimal_1beads/`

#### PhaseB1 検証: Torque あり，minimal_basal_stub (0.05s)

```bash
uv run python -m scripts.01_simulate_swimming \
  --duration-s 0.05 \
  --no-render-flagella \
  flagella.n_flagella=1 \
  flagella.stub_mode=minimal_basal_stub \
  motor.torque_Nm=4.0e-21 \
  motor.local_hook_scale=8.0 \
  output.base_dir=outputs/verify_phaseb1_torque4e21_hook8
```

#### PhaseB1 検証: Full flagella (0.05s)

```bash
uv run python -m scripts.01_simulate_swimming \
  --duration-s 0.05 \
  --no-render-flagella \
  flagella.n_flagella=1 \
  flagella.stub_mode=full_flagella \
  motor.torque_Nm=4.0e-21 \
  motor.local_hook_scale=8.0 \
  output.base_dir=outputs/verify_phaseb1_torque4e21_full_flagella
```

#### Sweep 再現テスト（単一条件）

```bash
uv run python -m scripts.01_simulate_swimming \
  --duration-s 0.05 \
  --no-render-flagella \
  flagella.n_flagella=1 \
  flagella.stub_mode=minimal_basal_stub \
  motor.torque_Nm=3.0e-21 \
  motor.local_hook_scale=4.0 \
  output.base_dir=outputs/verify_sweep_torque3e21_hook4
```

#### マルチべん毛テスト（PhaseC 向け）

```bash
uv run python -m scripts.01_simulate_swimming \
  --duration-s 0.05 \
  --no-render-flagella \
  flagella.n_flagella=3 \
  flagella.stub_mode=full_flagella \
  motor.torque_Nm=2.5e-21 \
  motor.local_hook_scale=8.0 \
  output.base_dir=outputs/verify_multi_n3_torque2.5e21
```

### 出力ファイル

各実行は `output.base_dir/{タイムスタンプ}/` の下に以下を生成：

```
step_summary.csv          # 全ステップの診断・判定結果（source of truth）
body_diagnostics.csv      # Body shape 詳細メトリクス
manifest.json             # 実行条件・メタデータ
run.log                   # stdout/stderr ログ
swim_state_projection.csv # 投影後の状態（参考）
swim_state_full.npy       # 全状態（オプション）
3d/swim_{frame}.npy       # 3D フレーム（--render-flagella 時）
```

---

## セクション 2: Sweep 実行（run_motor_scale_sweep.py）

### 基本構文

```bash
uv run python -m scripts.run_motor_scale_sweep \
  [--param-name PARAM] \
  [--param-values VAL1 VAL2 ...] \
  [--sweep-param SWEEP_PARAM] \
  [--sweep-values VAL1 VAL2 ...] \
  [--duration-s DURATION] \
  [--n-flagella N] \
  [--stub-mode MODE] \
  [--output-dir DIR]
```

### パラメータ説明

| オプション | 型 | デフォルト | 説明 |
|-----------|-----|----------|------|
| `--param-name` | str | `motor.torque_Nm` | sweep の行軸パラメータ |
| `--param-values` | list | `[1.0e-21, ..., 4.0e-21]` | 行軸の値リスト |
| `--sweep-param` | str | `motor.local_hook_scale` | sweep の列軸パラメータ |
| `--sweep-values` | list | `[1.0, 2.0, 4.0, 8.0, 16.0]` | 列軸の値リスト |
| `--duration-s` | float | 0.05 | 各実行の時間 |
| `--n-flagella` | int | 1 | べん毛数 |
| `--stub-mode` | str | `minimal_basal_stub` | stub 構成 |
| `--output-dir` | str | `outputs/` | 出力ディレクトリ |

### よく使う Sweep 例

#### local_hook_scale sweep（2軸：torque vs hook_scale）

```bash
uv run python -m scripts.run_motor_scale_sweep \
  --param-name motor.torque_Nm \
  --param-values 1.0e-21 2.0e-21 3.0e-21 4.0e-21 \
  --sweep-param motor.local_hook_scale \
  --sweep-values 1.0 2.0 4.0 8.0 16.0 \
  --duration-s 0.05 \
  --n-flagella 1 \
  --stub-mode minimal_basal_stub \
  --output-dir outputs/sweep_hook_scale_2026_04_16
```

出力: `outputs/sweep_hook_scale_2026_04_16/{タイムスタンプ}/`

#### local_spring_scale sweep（2軸：torque vs spring_scale）

```bash
uv run python -m scripts.run_motor_scale_sweep \
  --param-name motor.torque_Nm \
  --param-values 1.0e-21 2.0e-21 3.0e-21 4.0e-21 \
  --sweep-param motor.local_spring_scale \
  --sweep-values 0.5 0.8 1.0 1.5 2.0 \
  --duration-s 0.05 \
  --n-flagella 1 \
  --stub-mode minimal_basal_stub \
  --output-dir outputs/sweep_spring_scale_2026_04_16
```

#### local_bend_scale sweep（2軸：torque vs bend_scale）

```bash
uv run python -m scripts.run_motor_scale_sweep \
  --param-name motor.torque_Nm \
  --param-values 1.0e-21 2.0e-21 3.0e-21 4.0e-21 \
  --sweep-param motor.local_bend_scale \
  --sweep-values 0.5 1.0 2.0 4.0 \
  --duration-s 0.05 \
  --n-flagella 1 \
  --stub-mode minimal_basal_stub \
  --output-dir outputs/sweep_bend_scale_2026_04_16
```

#### local_torsion_scale sweep（2軸：torque vs torsion_scale）

```bash
uv run python -m scripts.run_motor_scale_sweep \
  --param-name motor.torque_Nm \
  --param-values 1.0e-21 2.0e-21 3.0e-21 4.0e-21 \
  --sweep-param motor.local_torsion_scale \
  --sweep-values 0.0 0.5 1.0 2.0 4.0 \
  --duration-s 0.05 \
  --n-flagella 1 \
  --stub-mode minimal_basal_stub \
  --output-dir outputs/sweep_torsion_scale_2026_04_16
```

### 出力ファイル

各 sweep 実行は `--output-dir/{タイムスタンプ}/` の下に以下を生成：

```
{sweep_name}_sweep_summary.csv          # 集約結果（body/hook/flagella pass/fail + metrics）
{sweep_name}_combined_pass_fail_heatmap.png    # メイン出力（1×3 subplot）
{sweep_name}_category_heatmap.png       # 診断軸（first_fail_category）
{sweep_name}_flagella_heatmap.png       # flagella pass/fail（データ有時のみ）
```

---

## セクション 3: ビジュアライゼーション（plot_motor_scale_collapse_heatmap.py）

### 基本構文

```bash
uv run python -m scripts.plot_motor_scale_collapse_heatmap \
  SWEEP_SUMMARY_CSV \
  [--output-dir DIR]
```

### パラメータ説明

| 引数 | 型 | 説明 |
|-----|-----|------|
| `SWEEP_SUMMARY_CSV` | path | sweep が生成した `*_sweep_summary.csv` |
| `--output-dir` | str | 出力先（デフォルト：CSV の親ディレクトリ） |

### よく使う可視化例

#### 直前の sweep 結果を可視化

```bash
uv run python -m scripts.plot_motor_scale_collapse_heatmap \
  outputs/sweep_hook_scale_2026_04_16/{TIMESTAMP}/hook_scale_sweep_summary.csv
```

出力: `outputs/sweep_hook_scale_2026_04_16/{TIMESTAMP}/` に heatmap 群を生成

#### 既存の sweep 結果を再度可視化

```bash
uv run python -m scripts.plot_motor_scale_collapse_heatmap \
  outputs/2026-04-15/234606/hook_scale_sweep_summary.csv \
  --output-dir outputs/analysis_2026_04_16/
```

出力: `outputs/analysis_2026_04_16/` に同じ heatmap を生成

---

## セクション 4: 等価性検証ワークフロー

### 目的

sweep 側の特定条件が単一実行側でも同じ動作をするか（step_summary.csv が同一か）を確認。

### 検証チェックリスト

以下を確認してから「等価である」と判定：

- [ ] **Configuration match**
  - [ ] `n_flagella` が一致
  - [ ] `stub_mode` が一致
  - [ ] `motor.torque_Nm` が一致（数値精度同一）
  - [ ] `motor.local_hook_scale` が一致（他の scale も同様）
  - [ ] `duration_s` が一致
  - [ ] `dt_s` が指定されている場合、両方に一致

- [ ] **Output comparison**
  - [ ] `step_summary.csv` の行数（step 数）が同じ
  - [ ] `first_fail_category` の値が同じ
  - [ ] `first_fail_step` が同じ

- [ ] **Diagnostic columns match**
  - [ ] `finite_pass`, `shape_pass_nonbody`, `body_shape_pass` の値系列が同じ
  - [ ] 回転追跡列（`flag_phase_rate_hz` 等）の NaN/Inf パターンが同じ

- [ ] **Body/hook metrics comparison**
  - [ ] `body_spring_stretch`, `body_bend_error` などの late-time 値がおおむね同じ
  - [ ] `hook_len_mean_over_b`, `local_attach_first_rel_err` が同じ傾向

### 実行例：Sweep から単一実行への検証

#### Step 1: Sweep で目視した条件を特定

例）`outputs/sweep_hook_scale_2026_04_16/{timestamp}/hook_scale_sweep_summary.csv` を確認し、
「torque=3.0e-21, local_hook_scale=4.0 で first_fail_category='none' が出た」と判定

#### Step 2: 同一条件で単一実行

```bash
uv run python -m scripts.01_simulate_swimming \
  --duration-s 0.05 \
  --no-render-flagella \
  flagella.n_flagella=1 \
  flagella.stub_mode=minimal_basal_stub \
  motor.torque_Nm=3.0e-21 \
  motor.local_hook_scale=4.0 \
  output.base_dir=outputs/verify_equiv_torque3e21_hook4
```

#### Step 3: 出力を比較

```bash
# sweep から該当行を抽出（例）
grep "torque.*3.0e-21.*hook.*4.0" outputs/sweep_hook_scale_2026_04_16/{timestamp}/hook_scale_sweep_summary.csv

# 単一実行の結果
cat outputs/verify_equiv_torque3e21_hook4/{timestamp}/step_summary.csv | tail -1

# 比較：first_fail_category と step 数が同じか確認
```

#### Step 4: 検証結果を記録

一致した場合:

```
✓ Verified equivalence: torque=3.0e-21, hook_scale=4.0
  - sweep step count: 500
  - single run step count: 500
  - first_fail_category: none (both)
  - verified datetime: 2026-04-16 10:30:00 UTC
  - run: outputs/verify_equiv_torque3e21_hook4
```

一致しない場合:

```
✗ Equivalence FAILED: torque=3.0e-21, hook_scale=4.0
  - sweep step count: 500
  - single run step count: 468 ← MISMATCH
  - first_fail_category: none (sweep) vs hook_break (single)
  - likely cause: different dt_s or random seed
  - TO DEBUG: compare dt_s and check rng state
```

---

## 参考リンク

- [prompts/phase2/0037.md](../../prompts/phase2/0037.md) - 詳細な実装方針・等価性ルール
- [docs/phase2/sim_diagnostics.md](sim_diagnostics.md) - 出力列の説明
- [docs/phase2/2026-04-14_motor_scale_sweep_report.md](2026-04-14_motor_scale_sweep_report.md) - sweep 実行履歴・結果
