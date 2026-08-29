# Issue #215 analysis implementation log

- 5秒campaignを再実行せず、compact stagingとaxis-angle comparisonの再現可能なCLIを実装した。
- cs10 NASで保存済みarchiveからmotion_featuresを生成し、raw archiveを除く904 MBのcompact artifactをlocal referenceへ転送した。
- 36 conditionを確認し、strict PASSは31件。non-PASS 5件はすべてn=4のshape_nonbody failureだった。
- 0--2秒はn=1--3で2秒referenceと一致し、2--5秒ではbody--flagella axis angleの単調収束は確認しなかった。
