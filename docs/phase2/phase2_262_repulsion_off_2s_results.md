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

## べん毛のまとまりと一時的な長軸方向の停滞

ユーザーの追加質問に対し、両armの各`diagnostic_samples.csv`から後半1.5–2.0 sの50時点を同じ方法で集計した。2本以上の12配置では、OFF側のべん毛らせん軸間の最大角の中央値がONより大きい配置は7、小さい配置は5。べん毛束の平均半径の中央値が大きい配置は8、小さい配置は4。たとえば`nf03__slots012`は最大軸間角24.1→44.4°、平均半径0.51→0.95 µm、`nf04__slots0124`は22.8→59.4°、0.86→0.96 µmと、OFF側で幾何学的なまとまりが弱い。一方`nf04__slots0134`は52.5→29.0°、0.77→0.75 µm、`nf06__slots012345`は65.1→47.0°、0.95→0.86 µmで逆方向である。**排除体積OFFが全配置で束化を悪化させたとはいえない。** これらは幾何学的なproxyであり、接触・絡まりを含む位相的な束化の判定ではない。各配置は1 seedだけである。

`nf03__slots024`の0.95–1.05 sでは、初期方向からの角度plotはほぼ平坦だが、3D長軸方向の両端差はOFFで1.20°、ONで1.31°。同じarchiveの30菌体beadを連続時点ごとに剛体最小二乗回転で合わせ、長軸方向へ射影した菌体材料点の軸まわり回転はOFFで36.80°、ONで36.05°だった。この時間帯は菌体の**自転が続いている**。

ただし、ユーザーが指摘した**16τ = 0.64 s付近では別の現象**がある。30菌体beadの材料点回転を1 ms archiveから積分すると、`nf03__slots024`の長軸まわり正味回転は以下の通り。符号は同じ向きの長軸を基準とし、負値は短時間の逆回転を表す。

| 時間帯 | ON | OFF |
|---|---:|---:|
| 0.40–0.60 s | +66.76° | +66.82° |
| 0.64–0.80 s | −2.81° | +0.96° |
| 0.80–1.00 s | +62.56° | +64.56° |

0.64–0.80 sは正味回転だけが相殺されたわけではない。同区間の絶対回転量の和もON 7.58°、OFF 7.78°と小さく、概ね停止して一部が微小に逆回転した。個別の5断面リングで独立に計算しても同区間の正味回転はON −5.39～−0.15°、OFF −1.26～+2.77°であり、30 bead剛体近似だけの見かけではない。剛体fitの最大RMS残差はこの区間でON 0.00566 µm、OFF 0.00751 µm。これらは変形可能な菌体の**最良剛体近似**であって厳密な剛体角速度ではない。

13配置の同時刻0.64–0.80 sもON/OFFの双方で計算した。`nf03__slots024`以外にも減速する配置はあるが、`nf01__slots0`は両armとも約+30.63°、`nf03__slots012`はON +38.84°、OFF +40.38°で回転が続く。全配置に共通する一斉停止とはいえない。`nf03__slots024`ではmotor torqueは両armで`2.5e-20 N m`、switching OFFで、strict QCもPASS。排除体積設定の違いが停止の必須条件ではないが、力学的原因はこの比較だけでは特定できない。初期軸からの角度という単一スカラーは軸方向の方位角変化を隠し得るため、長軸の向きと材料点の自転を別々に読む。

追加のarchive-only図・集計CSV・method manifestは`outputs/2026-10-19/023732/issue262_followup_analysis/`に保存した。菌体材料点の全13配置時系列、区間集計、16τ図、5断面リングの照合、再解析スクリプトは同ディレクトリの`body_roll/`に置いた。これは受入閾値の追加や再simulationではない。

## レビュー入口と限界

- ローカルレビュー：`outputs/2026-10-19/023732/issue262_2s_review.html`。13配置のQC、3D方向成分・角度plot、0.1秒窓、固定camera 3D/2D replayを開ける。
- 数値結果：`outputs/2026-10-19/023732/issue262_axis_diagnostic/`内の`axis_timeseries.csv`、`axis_windows.csv`、`paired_axis_windows.csv`、`axis_summary.csv`、plot、manifest。
- 共通評価：`outputs/2026-10-19/023732/issue262_2s_evaluation/`内の`summary.csv`、`window_qc.csv`、初期geometry、replay、manifest。

長軸方向の変化には連続的な旋回や緩和も含まれ、周期的wobbleとの区別はこの分類だけではできない。旧有限差分torsionとlegacy frame、束化・軸整列の未解決点も残る。OFF単独結果から原因を断定せず、この結果をcanonical modelやdatasetの採択根拠にしない。物理解釈を含むmergeはユーザー判断を待つ。
