# #262補足：16τ付近の菌体自転鈍化・原因監査

## 範囲と結論の強さ

#255の排除体積ON対照と#262の全OFF、いずれも`nf03__slots024`の完成済み2秒archiveを再解析した。新規simulationは行っていない。1 msごとの30菌体beadで長軸まわり材料点回転を測り、10 ms間隔の既存診断と、保存姿勢における力・Rotne–Prager–Yamakawa mobilityの再評価を照合した。両armは旧有限差分torsion、legacy attachment frame、局所one-ring全ベクトル反作用、motor switching OFF、Brownian OFF、同一seedである。

**最有力の実装上の機序**は、べん毛根元の最初のsegmentがべん毛の主軸とほぼ平行になり、segment force coupleの有効な腕が縮むこと。一定の軸方向指令トルクから、過大な横方向トルクと力が生じる。その全ベクトル反作用が菌体の局所one-ringに渡り、幾何形状・弾性力・流体力学的相互作用との合成で長軸自転が鈍化する。時間的・数値的な整合性は強いが、修正モデルの**軌道を実行していないため因果確定ではない**。

## 時系列で確認した事実

| `nf03__slots024`の量 | 0.40–0.60 s | 0.64–0.80 s | 0.86–0.95 s |
|---|---:|---:|---:|
| ON菌体自転、40 ms平均の中央値 | +0.91 Hz | −0.063 Hz | +1.11 Hz |
| OFF菌体自転、同上 | +0.91 Hz | −0.015 Hz | +1.12 Hz |
| ON根元segmentの最小横方向lever arm中央値 | 0.351 µm | 0.0073 µm | 0.283 µm |
| OFF同上 | 0.351 µm | 0.0076 µm | 0.310 µm |
| ON適用motor torqueの1べん毛あたり平均norm／指令値 | 1.06 | 9.26 | 1.01 |
| OFF同上 | 1.06 | 9.70 | 1.01 |
| OFF最大横方向motor torque／指令値の中央値 | 0.37 | 11.37 | 0.13 |

自転の正味回転量は0.64–0.80 sにON −2.81°、OFF +0.96°で、0.80–1.00 sにはON +62.56°、OFF +64.56°へ戻る。5つの菌体断面リングを別々に追跡しても鈍化が出る。根元のlever arm縮小、motor torqueの増幅、自転鈍化は両armのほぼ同時刻に現れ、0.85 s頃にlever armが回復すると自転も戻る。ONのこの時間帯には小さな実repulsion forceが一部時点であるが、OFF側の実repulsion forceはゼロで同じ現象が起きる。`flag_flag_repulsion_force_*`は有効/無効に関係なく距離から計算する仮想force指標なので、実際に積分されたforceには`F_repulsion_mean_flag`等を用いた。

## 力の式と増幅の経路

`compute_root_torque_segment_couples_forces`では、べん毛主軸の単位ベクトルを`a`、segmentの両端差を`r`、軸に垂直な成分を`r⊥`、そのsegmentの配分重みを`w`とすると、力対の片側は

```text
F = τ w (a × r) / |r⊥|²
|F| = |τ w| / |r⊥|
|T_transverse| / |T_axial| = |r_parallel| / |r⊥|
```

となる。軸方向トルクは`τw`に合うが、`r⊥ → 0`で力と横方向トルクが増大する。今回の最初のflagellar segmentは0.64–0.80 sに`r⊥`がおよそ0.007 µmとなり、根元segmentはべん毛主軸に約1°以内で沿う。配分重みは同時間帯も約0.15でゼロではない。OFFの最大motor駆動力normの中央値は0.40–0.60 sの約0.027 pNから0.64–0.80 sの約0.69 pNへ増えた。モデルはこの全ベクトルトルクに等しく逆向きの反作用を菌体attachment周囲のone-ringに与える。流体中の速度は各beadの力の単純な軸成分ではなく、RPY mobility行列と全力の積で決まる。

保存姿勢に対して、同じ決定論的local twist状態を再構成し、力とmobilityを**時間を進めず**評価した。再構成した根元orientationは記録値と最大`3.3e-10°`以内で一致した。0.67 sのモデル瞬時菌体軸自転はON −3.30 Hz、OFF −2.73 Hz。motorの横方向トルクだけを`minimum_norm`補正で除いた固定姿勢反実仮想では、それぞれ−0.33 Hz、+0.90 Hzとなる。一方0.79 sのOFFは+0.79 Hzから−2.23 Hzへ変わり、全時点で改善しない。補正後も根元の小lever armから生じる大きな**局所force自体は残る**。この比較は完全な再計算軌道の代わりにはならない。

## 代替説明と残る疑問

- 指令motor torqueは各べん毛`2.5e-20 N m`で一定。switching、ramp、BrownianはOFFであり、16τの予定イベントはない。失速をmotor停止とは解釈できない。
- 旧有限差分torsionと同じエネルギーの解析的勾配を保存姿勢で比較すると、選んだ0.55–0.86 sの時点で菌体自転への差は最大0.027 Hz。有限差分誤差がこの瞬間の主因とは見えない。ただし解析的torsionで全軌道を再実行した比較ではない。
- 0.64–0.80 sのOFFでは10 ms診断での最大1 step変位が約0.010 µmに達する。有限値・形状のstrict QCはPASSだが、時間刻み収束と力の特異性の関係は未検証。
- 横方向トルク増大だけで停止は決まらない。同時刻の`nf06__slots012345`はOFFで横方向トルク比の中央値12.6でも菌体が+65.4°回る。各motorのベクトル和、attachment配置、弾性反力、流体力学的相互作用が長軸まわり速度を決める。
- 根元segmentが揃って主軸に近づく**上流の原因**は未特定。hook/attachment frame、初期位相の同期、弾性形状、局所twist分配の関与を分離する必要がある。
- 1τ（0.04 s）の事前screenは16τ（0.64 s）の現象を原理的に捉えられない。

## 次に識別力のある検証

1. 現行モデルの保存姿勢に対して、最初のsegmentの`|r⊥|`、駆動力norm、横方向トルク比、body/flagella各材料点回転を時間分解能を揃えて常時出す。まず診断値とし、新しい合否閾値は置かない。
2. 根元segmentに依存して`1/|r⊥|`が発散しないmotor force coupleを、反作用の局所性・力/トルク保存を保つ**別候補**として設計する。既存の`minimum_norm`横トルク補正だけでは大きな局所forceを消せない。
3. その候補とhook/attachment frame変更を混ぜずに、少なくとも25τ（1秒）のshort stageで失速・回復区間を覆い、時間刻み・位相seed・配置を変えて比較する。long-durationへ進む前にstrict QCと固定camera replayを確認する。cs10での新規実行には所定の個別承認が必要。

再解析入口は`outputs/2026-10-19/023732/issue262_followup_analysis/cause_analysis/`。`nf03_slots024_rotation_cause_timeline.png`、`nf03_slots024_event_timeline.csv`、`motor_snapshot_decomposition.csv`、`motor_per_flag_snapshot.csv`、method manifest、再実行スクリプトを置いた。これは#262の排除体積OFF分類とは別のモデル安定性監査であり、canonical採択やmerge判断を変更しない。
