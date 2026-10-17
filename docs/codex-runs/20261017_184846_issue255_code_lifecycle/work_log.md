# PR限定物の棚卸しと共通化

- ユーザー指示：merge前に削除・引継ぎ・共通化を区別し、今後も適用する。
- C単独2sは未使用・見送りのためcampaign/job/testを削除。復元先792d0d8を記録。
- 実行済み・取消済み予約の契約は保全し、結果との参照を維持。
- previewは実座標・接続で重複判定。代表はcampaign順、countはmodelから取得し、固定slot/反作用arm/表示範囲を除去。画像とmanifestを整合。
- 本数混在時の束診断CSV列欠落／例外を修正し、入力順を逆にした場合も検証。
- full pytest1004 PASS、targeted・lint・diff PASS。画像はtest一時出力で目視確認。
- 共通workflowとskill routerにmerge前の棚卸しを追加。新しいsimulation/cs10操作なし。
- Cloud review・mergeは引き続き保留。
