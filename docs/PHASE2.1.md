# Phase2.1 サマリ（データ生成 MVP）

## 達成したこと
- 剛体菌体＋螺旋べん毛（デフォルト5本）での遊泳シミュレーション骨格を実装（推進＋カウンター回転＋ブラウン運動）。
- 出力物（5秒・50fps）:
  - `sim/trajectory.csv`
  - `render/swim3d.mp4`, `render/swim3d_final.png`（Matplotlib 3Dグリッド表示）
  - `render2d/projection.mp4`, `render2d/frames/*.png`
  - `run.log`, `manifest.json`
- べん毛を螺旋形状で描画（色分け＋F0..ラベル）。菌体姿勢に追随。
- CLI上書き対応: `--duration-s`, `--fps-out`, `--render-flagella` と任意の `key=value` overrides。log/manifestは上書き後の設定を記録。
- tracking/* は Phase2 では出力しない方針に更新。

## 現状の挙動
- 推進速度: `v_prop = swim_coeff * n_flagella * f_motor * body_axis` （swim_coeff=0.0025）。
- 逆回転: `omega = -spin_coeff * f_motor * body_axis`（spin_coeff=0.05）。
- べん毛のモータ自転（軸回り回転）は未実装。べん毛は菌体と剛体結合のみ。
- FOV: 256px × 0.203µm/px ≈ 52µm立方体。初期位置を原点へ平行移動して中央配置。座標の正規化なし。

## 残タスク（Phase2.2以降）
1. 物理係数（swim_coeff / spin_coeff）の根拠づけとconfig露出、コメント整備。
2. manifestに最終適用configスナップショットと生成ファイル一覧を自動記録。
3. べん毛の自転（モータ回転位相）を描画に反映するか検討。
4. 速度スケール妥当化テスト（絶対値を実測レンジに合わせる）。
5. はみ出し対策オプション（クリップ/リセンタリング/ミニマップ）を設計。
6. 3D見た目調整（グリッド間隔・カメラ角度・凡例表示タイミング）。
