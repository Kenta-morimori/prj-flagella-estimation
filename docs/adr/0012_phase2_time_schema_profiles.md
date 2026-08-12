# ADR 0012: Phase 2 time schema と profile別時間スケール

- status: accepted
- date: 2026-08-02
- scope: Phase 2 / Issue #165

## Context

従来の Phase 2 config は `time.duration_s` と `time.dt_star` を直接持ち、`tau_s` は実装上 `1.0 s` 固定だった。これは 2010 project profile の過去run互換には必要だが、2015 refined model では論文の `dt_star=1e-5` と `tau = eta b^3 / T` に基づく実時間換算を分離して扱う必要がある。

また `time.dt_s` は名前上は内部積分刻みに見えるが、現行運用では出力・記録間隔として使われる。今回のPRでは互換性を優先し、意味を変えず deprecated key として残す。

## Decision

正式な時間schemaは次とする。

```yaml
time:
  duration:
    value: 1.0
    unit: tau  # tau / s
  integration:
    dt_star: 1.0e-5
```

`time.duration_s`、`time.dt_star`、`time.dt_s`、CLI shorthand `--duration-s` は deprecated compatibility として受け付ける。新schemaとlegacy flat keyが同時に指定され、同じ意味の値が矛盾する場合は警告ではなく `ValueError` で拒否する。同値なら `time_schema_source=mixed_equivalent` として受理する。

profile別の時間スケールは次のように固定する。

| profile | `tau_s` policy |
| --- | --- |
| `conf/sim_swim_2010.yaml` | legacy fixed `tau_s=1.0` |
| `conf/sim_swim_2010_paper.yaml` | `eta * b_m^3 / abs(reference_torque_Nm)` |
| `conf/sim_swim_2015.yaml` | `eta * b_m^3 / abs(reference_torque_Nm)` |
| `conf/sim_swim_2015_paper.yaml` | `eta * b_m^3 / abs(reference_torque_Nm)` |

`motor.enabled` は実駆動ON/OFF、`motor.torque_Nm` は実際に加える符号付きtorque、`motor.reference_torque_Nm` は時間換算と無次元係数の参照torqueを表す。motor-offでも `reference_torque_Nm` は保持する。

`motor.torque_Nm=-1` sentinel は 2010 paper profile由来の legacy compatibility としてだけ `reference_torque_Nm = eta b^3`, `motor_torque_Nm = eta b^3`, `motor_enabled = true` に正規化する。profileなしの古いfixtureは parser互換のため維持するが、canonical 2010 project / 2015 profileでは拒否する。

## Step 時刻契約

内部step数は次で決める。

```text
total_steps = ceil(duration_tau / dt_star)
```

出力契約は次のとおり。

```text
initial state:
  t = 0

number of states:
  total_steps + 1

step_summary.csv rows:
  total_steps

final step_summary time:
  (total_steps - 1) * dt_star * tau_s

final state time:
  total_steps * dt_star * tau_s
```

`step_summary.csv` は各stepの開始時刻を記録する。`ceil` により final state時刻は指定durationをわずかに超える場合がある。この超過は内部積分契約の結果であり、manifestに `final_state_t_s` と `final_step_summary_t_s` を保存して追跡する。

## Manifest

user-facing logとmanifestの論文表記は、次を優先する。

```text
t       実時間 [s]
τ       時間尺度 [s]
Δt      実時間の積分刻み [s]
t/τ     無次元時刻 [-]
Δt/τ    無次元積分刻み [-]
```

machine-readable manifestでは `time.paper_notation` にASCII keyの `t`、`tau`、`delta_t`、`t_over_tau`、`delta_t_over_tau` を保存する。`t_star`、`dt_star`、`duration_tau`、`tau_s`、`dt_internal_s` は既存分析・runとの互換keyとして維持する。`dt`単独は、実時間刻みと無次元刻みを区別できないため、user-facing log・文書では使用しない。

single-run manifest と generic multi-run `run_manifest.json` には、正規化済み時間量として次を保存する。

```text
duration_value
duration_unit
duration_tau
duration_s
dt_star
tau_s
dt_internal_s
total_steps
final_state_t_s
final_step_summary_t_s
time_scale_policy
reference_torque_Nm
motor_enabled
motor_torque_Nm
time_schema_source
legacy_time_keys_used
```

`time_schema_source` は `canonical`, `legacy_duration_s`, `legacy_dt_star`, `cli_shorthand`, `mixed_equivalent` のいずれかとする。`legacy_time_keys_used` は旧campaignやreplayがどの互換経路を通ったかを追跡するために保存する。

## Boundaries

2015 profilesは引き続き `implementation_status: pending` とし、parse・正規化・manifest出力の対象に留める。30-bead geometry、2015 paper motor dynamics、stability評価はそれぞれ Issue #166 / #167 / #168 で扱う。

`time.dt_s` から `output.interval_s` または `time.output_interval_s` への名称移行は後続Issueで扱う。今回のPRでは既存互換を優先し、`time.dt_s` の意味を変更しない。

## Verification status

parser unit test、canonical profile test、CLI manifest test、generic multi-run manifest test、invalid schema test、simulation loop契約の既存testで確認する。長時間simulationとvisual review動画生成はIssue #165の完了条件に含めない。
