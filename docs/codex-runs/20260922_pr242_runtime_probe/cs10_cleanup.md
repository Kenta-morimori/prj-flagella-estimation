# PR #242 cs10 cleanup（2026-09-22 JST）

- queueをpause後、reservation 5を正式cancelした。最終stateは`cancelled`、exit code `-15`、通知試行済み1回。予約6は`queued`のまま。
- direct jobのPGID 3328（nf4/nf5）にSIGTERMを送った。reservation 5のPGID 45559（nf6を含む）も終了し、対象processが存在しないことを確認した。
- 3 childの`campaign_completion.json`はいずれも`status=running`で、`run_summary.json`、`performance.json`、compact checkpointは存在しなかった。削除前の大きさはnf6 28K、nf4 24K、nf5 28K。

以下の未完走childディレクトリだけを削除し、不在を確認した。NAS上の通常削除であり、復元は管理者側のバックアップが存在する場合に限られる。

1. `/net/fs01/volume1/work01/Ktakemori/prj-flagella-estimation/outputs/2026-09-08/044446/parallel/queue-00005-14c7fb53__ccff7345439f/children/006_nf06`
2. `/net/fs01/volume1/work01/Ktakemori/prj-flagella-estimation/outputs/2026-09-15/140131/parallel/issue184_2015_nf4_nf5_10tau__direct2/children/001_nf04`
3. `/net/fs01/volume1/work01/Ktakemori/prj-flagella-estimation/outputs/2026-09-15/140131/parallel/issue184_2015_nf4_nf5_10tau__direct2/children/002_nf05`

旧reservation 5の完走済みnf1–3、親jobのmanifest、停止済みdirect nf1/nf2 rootは保持した。旧childを新benchmarkへコピーしない。
