# Repo Rules（実装・運用ルール）

このファイルは「毎回共通の制約」を定義する。タスク固有の要求（DoD、作業手順、削除対象の列挙など）は `prompts/Phase*/` に書く。

## 1. scripts / src の責務分離（必須）
- `scripts/*.py` は CLI入口である。責務は以下に限定する：
  - config読み込み（YAML + key=value override）
  - 出力ディレクトリ作成
  - 実行コンテキスト初期化（ログ/manifest作成の統一）
  - `src` の関数呼び出し（例：`pipeline.run(config, out_dir)`）
- アルゴリズム/処理本体は必ず `src/` に置く（scriptsにロジックを増やさない）。

## 2. 実行コンテキスト（init_run）
- `init_run()` は「実行の器」を統一するための共通関数である：
  - `outputs/YYYY-MM-DD/HHMMSS/` の作成（JST）
  - `run.log` と `manifest.json` の作成
  - git clean 強制と commit id 記録（有効化する場合は一貫して適用）
- 各CLIは原則 `init_run()` を通す（例外を作らない）。

## 3. 出力規則（必須）
- 実行結果は必ず `outputs/YYYY-MM-DD/HHMMSS/` に保存する（JST）。
- `run.log` と `manifest.json` を必ず生成する。
- `manifest.json` には少なくとも以下を記録する：
  - configパス、override（あれば）、seed（あれば）
  - 入出力パス
  - 実行環境情報（git commit / python version等、記録可能な範囲）

## 4. パッケージ命名と依存（必須）
- Phase2実装のトップレベルパッケージ名は `sim_swim` とする。
- Phase3-4実装のトップレベルパッケージ名は `flagella_estimation` とする。
- importは上記のトップレベルから開始し、相対importの乱用は避ける。

## 5. 例外・ログ・品質（必須）
- 例外は握り潰さず、原因が分かるメッセージを付与して送出する。
- 処理の開始/終了、主要パラメータ、入出力パス、件数、所要時間をログに残す。
- 公開関数には日本語docstring（概要 + Args/Returns）を付与する。

## 6. テスト方針（恒久ルール）
- 「走ること」を担保する最小のスモークテストは維持する（主要CLIまたは主要pipeline）。
- `init_run()` と出力規則（run.log/manifest生成）はテストまたは簡易検証で担保する。
- それ以上の詳細テスト（物理妥当性、モデル精度など）は Phase別プロンプトでDoDとして定義する。