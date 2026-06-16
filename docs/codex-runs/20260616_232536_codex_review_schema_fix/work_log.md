# Codex review schema fix

## 目的

PRコメント `@codex review` で起動する `codex-review` workflow が、OpenAI structured output schemaの検証で失敗する問題を修正する。

## 原因

`inline_comments.items.properties` には `suggested_code` が定義されていたが、`additionalProperties: false` のobject schemaで `required` に含まれていなかった。
Codex/OpenAI API側では、structured output用schemaとして各objectの `required` が `properties` の全キーを含む必要があり、`Missing 'suggested_code'` で失敗していた。

## 変更

- `.github/workflows/codex-review.yml` の `inline_comments.items.required` に `suggested_code` を追加した。
- `suggested_code` がない場合は空文字を返すよう、Codex review promptに明記した。
- ADR 0005にschema上の必須キーと空文字運用を追記した。

## 確認

- YAML構文確認: `ruby -e 'require "yaml"; YAML.load_file(".github/workflows/codex-review.yml"); puts "YAML OK"'`
- output-schema JSON確認: `SCHEMA OK`
- git diff確認: `git diff --check`
