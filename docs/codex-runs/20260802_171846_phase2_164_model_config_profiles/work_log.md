# Phase 2 Issue #164 work log

## Scope

2010年・2015年のproject/paper model configを4 profileへ再編し、profile identityとsource config pathをmanifestから追跡可能にした。Issue #165の時間scale契約、#166のrefined geometry、#167/#168の2015 dynamics・評価は実装していない。

## Implementation

- 旧 `conf/sim_swim.yaml` を、物理・数値設定を維持した `conf/sim_swim_2010.yaml` へrenameした。
- 2010 paper参照profileと2015 refined project/paper pending profileを追加した。
- typed `model_profile`、nominal bead数整合、enum検証、provenance override拒否を追加した。
- `implementation_status=pending` の実simulationをoutput作成前に拒否するようにした。
- single-runとgeneric multi-run、および既存sweep manifestへ `model_profile` / `source_config_path` を追加した。
- CLI、campaign config、sweep、replay fallback、active docs/testsのdefault pathを `conf/sim_swim_2010.yaml` へ更新した。
- ADR 0010へprofile構成、paper/project差分、後続Issueとの境界を記録した。

## Verification

- focused tests: 164 passed
- light tests: 246 passed, 240 deselected
- full pytest: 493 passed
- Ruff format/check: PASS
- all YAML parse: PASS
- `git diff --check`: PASS
- build-only check: 2010 project 48 beads、2010 paper 60 beads、2015 project/paperはpending error
- 2010 project CLI 1-step smoke: PASS、`/private/tmp/phase2_164_profile_smoke/manifest.json` にprofile/source pathを確認
- executable legacy path search: source、scripts、tests、conf、active READMEで0件

長時間simulationと動画生成は実行していない。config/profile/manifest契約の変更であり、ユーザーvisual reviewは不要と判断した。

## Codex review follow-up

Codex Cloud reviewは、YAMLの非整数profile値が`int()`で切り捨てられる点を指摘した。`year`と4つのbead/count fieldについて、bool・float・文字列を拒否する厳密な整数validationへ変更し、6ケースの回帰testを追加した。修正後のfocused/full testとRuffを再実行した。
