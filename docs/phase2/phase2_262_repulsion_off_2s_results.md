# Issue #262：排除体積全OFF・2秒診断

## 実行と入力の照合

cs10 queue #28は固定commit `a444a1bd45422d9b5b911df44ee61b397a30f809`、`conf/phase2_parallel/issue262_repulsion_off/hex_2s_job.yaml`で2026-10-17 22:20頃に開始し、2026-10-19 02:33 JSTに最後の配置が完走した。13配置すべてが`50τ = 2.0 s`、500,000 stepを連続実行し、1 ms間隔のarchiveを保存した。`campaign_completion.json`は`completed_condition_count=13`、queueは`succeeded`、各`run_summary.json`は`execution.status=completed`、`row_count=expected_total_steps=500000`、連続step、online finite/shape PASSである。

ローカルの`outputs/2026-10-19/023732/issue262_hex_repulsion_off_2s/`にportable artifact 122件を同期し、`sync_verification.json`でcs10とのSHA-256全件一致を確認した。別途13件の子`run_manifest.json`も`issue262_condition_manifests/`に同期し、`condition_manifest_verification.json`でSHA-256を確認した。13件ともbody–flagella、flagella–flagellaのcandidateは存在し、`body_flagella_active=flagella_flagella_active=active_total=0`、両enable flagはfalseである。

共通`development_evaluation`を#262元設定で実行し、strict finite/body/hook length/flag/motorは**13/13 PASS**。初期geometryと固定camera 3D/2D replayを生成した。replayは本数別の3D/2D各6本、いずれも41 frameで全12本をdecode確認した。対照の#255局所反作用ON archiveは元の`conf/phase2_multi_run/2010_hex_project_local_reaction_2s_issue255.yaml`で別評価済み（13/13 strict PASS、13/13完走、各500,000 step、122/122 SHA一致）。#262と対照を同じ評価入力へ混ぜていない。

`body_axis_diagnostic`は13 ID、base config、条件別overrides、0.1秒窓を照合し、両設定の物理的差が`potentials.spring_spring_repulsion.flagella_flagella_enabled=false`だけであることを確認してから同時刻比較を作成した。対照との比較はこの旧有限差分torsionモデル内に限る。

## 菌体長軸の3D方向変動

`axis_timeseries.csv`は各配置の符号付き3D単位方向`axis_x/y/z`、初期方向からの角度、連続archive間の角度を含む。各配置のcompleted archiveは2002 statesを持ち、`axis_windows.csv`は0.1秒ずつ20窓を記録する。下表は方向の時間変化が最も大きいOFF側0.1秒窓を示す。**「揺れあり」は2秒内に3D方向変動が観測されたという記述的分類であり、周期性や原因を判定する新閾値ではない。** 全13配置で変動が観測され、定常0、判定不能0と分類する。

| 配置 | 分類 | 最大初期角 | OFF最大0.1秒端点変動 | 該当時間帯 | ON最大0.1秒端点変動 |
|---|---|---:|---:|---|---:|
| `nf01__slots0` | 揺れあり | 37.9° | 15.9° | 0.7–0.8 s | 15.9° |
| `nf02__slots01` | 揺れあり | 85.1° | 21.6° | 0.7–0.8 s | 19.8° |
| `nf02__slots02` | 揺れあり | 60.2° | 17.9° | 1.5–1.6 s | 13.5° |
| `nf02__slots03` | 揺れあり | 25.8° | 11.9° | 1.2–1.3 s | 16.7° |
| `nf03__slots012` | 揺れあり | 107.9° | 26.9° | 0.8–0.9 s | 26.0° |
| `nf03__slots013` | 揺れあり | 69.5° | 19.9° | 1.3–1.4 s | 18.8° |
| `nf03__slots014` | 揺れあり | 60.2° | 27.1° | 0.7–0.8 s | 22.1° |
| `nf03__slots024` | 揺れあり | 27.1° | 16.2° | 1.6–1.7 s | 17.1° |
| `nf04__slots0123` | 揺れあり | 91.6° | 30.1° | 0.3–0.4 s | 31.1° |
| `nf04__slots0124` | 揺れあり | 83.0° | 30.6° | 1.1–1.2 s | 33.6° |
| `nf04__slots0134` | 揺れあり | 35.3° | 20.0° | 1.9–2.0 s | 29.8° |
| `nf05__slots01234` | 揺れあり | 80.1° | 41.6° | 0.3–0.4 s | 39.2° |
| `nf06__slots012345` | 揺れあり | 84.1° | 35.8° | 1.6–1.7 s | 43.4° |

対照の最大値は別の窓で生じる場合があり、表の最大値同士は同時刻差ではない。同時刻の比較は`paired_axis_windows.csv`の260窓を使う。ONより変動が小さい窓も大きい窓もあり、排除体積全OFFで方向変動が一様に消えたとはいえない。1本配置のON/OFF時系列は一致し、べん毛間pairを切る操作の内部対照になっている。

## レビュー入口と限界

- ローカルレビュー：`outputs/2026-10-19/023732/issue262_2s_review.html`。13配置のQC、3D方向成分・角度plot、0.1秒窓、固定camera 3D/2D replayを開ける。
- 数値結果：`outputs/2026-10-19/023732/issue262_axis_diagnostic/`内の`axis_timeseries.csv`、`axis_windows.csv`、`paired_axis_windows.csv`、`axis_summary.csv`、plot、manifest。
- 共通評価：`outputs/2026-10-19/023732/issue262_2s_evaluation/`内の`summary.csv`、`window_qc.csv`、初期geometry、replay、manifest。

長軸方向の変化には連続的な旋回や緩和も含まれ、周期的wobbleとの区別はこの分類だけではできない。旧有限差分torsionとlegacy frame、束化・軸整列の未解決点も残る。OFF単独結果から原因を断定せず、この結果をcanonical modelやdatasetの採択根拠にしない。物理解釈を含むmergeはユーザー判断を待つ。
