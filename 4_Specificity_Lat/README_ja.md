# 4_Specificity_Lat：緯度ブロック化randomization

`3_Specificity_SpatialBlocks` と同じsampling-unit label shuffleを、continentではなく符号付き緯度帯で制約します。北緯と南緯は分け、デフォルトでは20度幅です。

同じ `sample_id` / `sample_ids` を共有するsampling unit間では植物ラベルを交換しません。

制約付きshuffleは、sample_id非互換ペアを事前計算し、通常shuffle後に違反割り当てだけを修復する高速実装です。

```bash
cd "[CURRENT DIRECTORY]"
MPLCONFIGDIR=/tmp/mplconfig python3 latitude_blocks_sh_workflow.py \
  --globalfungi-results ../1_GlobalFungi/results_helotiales \
  --output-dir outputs \
  --plant-rank family \
  --fungi-selection rank --fungi-rank-max 150 \
  --plant-selection rank --plant-rank-max 30 \
  --latitude-block-width 20 \
  --min-block-unique-occurrences 100 \
  --n-randomizations 10000 \
  --jobs 0 \
  --offline
```

真菌属版の例です。

```bash
cd "[CURRENT DIRECTORY]"
MPLCONFIGDIR=/tmp/mplconfig python3 latitude_blocks_genus_workflow.py \
  --globalfungi-results ../1_GlobalFungi/results_helotiales \
  --output-dir outputs \
  --plant-rank family \
  --fungi-selection rank --fungi-rank-max 150 \
  --plant-selection rank --plant-rank-max 30 \
  --latitude-block-width 20 \
  --min-block-unique-occurrences 100 \
  --n-randomizations 10000 \
  --jobs 0 \
  --offline
```

このステップの出力run名を、最後の `5_Phylogenetic_Signal --specificity-run-name` に渡します。
