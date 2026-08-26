# 2_Specificity：blockなしrandomization

`1_GlobalFungi/results_helotiales` の medium-or-higher 出力を読み込み、run内では `helotiales_root_occurrences_with_hosts.tsv` という汎用名にコピーして解析します。

ランダマイゼーションではsampling unitの植物ラベルおよびcontinentラベルをシャッフルします。同じ `sample_id` / `sample_ids` を共有するsampling unit間ではラベル交換を禁止しています。

制約付きshuffleは、sample_id非互換ペアを事前計算し、通常shuffle後に違反割り当てだけを修復する高速実装です。

```bash
cd "[CURRENT DIRECTORY]"
MPLCONFIGDIR=/tmp/mplconfig python3 specificity_sh_workflow.py \
  --globalfungi-results ../1_GlobalFungi/results_helotiales \
  --output-dir outputs \
  --plant-rank family \
  --fungi-selection rank --fungi-rank-max 150 \
  --plant-selection rank --plant-rank-max 30 \
  --n-randomizations 10000 \
  --jobs 0 \
  --offline
```

真菌属版の例です。

```bash
cd "[CURRENT DIRECTORY]"
MPLCONFIGDIR=/tmp/mplconfig python3 specificity_genus_workflow.py \
  --globalfungi-results ../1_GlobalFungi/results_helotiales \
  --output-dir outputs \
  --plant-rank family \
  --fungi-selection rank --fungi-rank-max 150 \
  --plant-selection rank --plant-rank-max 30 \
  --n-randomizations 10000 \
  --jobs 0 \
  --offline
```

`--jobs 0` は利用可能CPU数を使います。
