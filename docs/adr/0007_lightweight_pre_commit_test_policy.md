# ADR 0007: Lightweight pre-commit test policy

## Status

Accepted

## Context

従来の `.githooks/pre-commit` は commit ごとに以下を実行していた。

* `uv run ruff format --check .`
* `uv run ruff check .`
* `uv run pytest -q`

この設定は CI と同等の safety を局所的に得やすい一方で，Codex CLI のように小さな修正を高頻度で commit する運用では待ち時間とログ確認が増えやすい。Phase 2 には simulation，render，sweep，dataset，CLI 系の比較的重い pytest も含まれており，毎 commit で full pytest を回すと開発速度を不必要に落とす。

一方で full pytest 自体は，モデル変更や出力仕様変更の検証として引き続き必要であり，CI でも維持されている。

## Decision

pre-commit hook の既定を lightweight checks に変更する。

既定:

* `ruff format --check .`
* `ruff check .`
* `pytest -q -m light`

full pytest は削除せず，以下で実行する。

* モデル・physics・geometry・simulation core の変更時
* output format / manifest / CSV schema の変更時
* dataset 生成仕様の変更時
* PR 前
* CI

初回 `light` 対象は，短時間・deterministic・subprocess/render/sweep を含まない unit test に限定する。include 方式を採用し，`pytest.mark.light` を付与した test のみを pre-commit の既定対象にする。

また，必要に応じて hook 実行時に `FULL_TEST=1` を指定すれば full pytest を明示実行できるようにする。

## Consequences

利点:

* 小さな commit の待ち時間を減らせる
* Codex CLI の不要なトークン消費を抑えられる
* full pytest は PR 前と CI で維持される

注意点:

* pre-commit を通っただけでは full regression の代替にならない
* `light` marker の対象選定は保守が必要
* 高リスク変更時に full pytest を省略しない運用徹底が必要
