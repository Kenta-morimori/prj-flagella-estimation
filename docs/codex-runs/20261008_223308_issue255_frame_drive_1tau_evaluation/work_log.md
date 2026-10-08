# Issue #255 queue #22 1tau評価

cs10読み取り接続でqueue #22成功・#21 cancelledを確認。3条件が完走。portable artifactをsync_reference_from_cs10.py parallel-campaignで同期し32ファイルSHA-256一致を記録。run_summaryを先に読み共通evaluatorで3/3 PASSと固定camera6動画を生成。動画は41フレームを確認し3D終端・2D配置・時系列を目視した。旧#20の同じ配置を0～0.04秒へ揃えて物理設定・初期位置・source commitを照合して比較した。全step motor/frame釣り合いの修正を確認したが密な束化は未確認、nf5の遊泳方向もずれている。結果文書へ1s候補3条件・見積り約9時間26分を記録。1s予約・開始、#21の追加操作、canonical採択はしない。source code変更なしでruntime/QC/replayの結果検証を優先しfull pytestは繰り返さない。
