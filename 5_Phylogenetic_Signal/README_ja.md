# 5_Phylogenetic_Signal：系統シグナル解析

このステップは最後に実行します。デフォルトでは `4_Specificity_Lat` のrunフォルダから、宿主特異性の統合表、2DP表、FASTAを自動取得します。

`3_Specificity_SpatialBlocks` の結果を使う場合は、`--specificity-output-root ../3_Specificity_SpatialBlocks/outputs` を明示してください。

## 4_Specificity_Lat の結果を使う例

```bash
cd "[CURRENT DIRECTORY]"
python3 phylogeny_workflow.py \
  --data-dir data \
  --specificity-run-name sh150_family30_10000_lat20_minblock100 \
  --specificity-fasta-subset its1_and_its2 \
  --target-taxon Helotiales \
  --trait-column host_z_standardized_dprime_sampling_unit_label_shuffle \
  --minimum-length 400 \
  --outgroup-minimum-length 460 \
  --bootstrap 500 \
  --signal-randomizations 10000 \
  --jobs 0 \
  --iqtree-thread-mode auto \
  --figure-font-scale 1.2
```

## 3_Specificity_SpatialBlocks の結果を使う例

```bash
cd "[CURRENT DIRECTORY]"
python3 phylogeny_workflow.py \
  --data-dir data \
  --specificity-output-root ../3_Specificity_SpatialBlocks/outputs \
  --specificity-run-name sh150_family30_10000 \
  --specificity-fasta-subset its1_and_its2 \
  --target-taxon Helotiales \
  --trait-column host_z_standardized_dprime_sampling_unit_label_shuffle \
  --minimum-length 431 \
  --outgroup-minimum-length 460 \
  --bootstrap 1000 \
  --signal-randomizations 10000 \
  --jobs 0 \
  --iqtree-thread-mode auto \
  --figure-font-scale 1.2
```

`--minimum-length` は入力ターゲット配列の下限長、`--outgroup-minimum-length` はUNITEから探索する外群候補配列の下限長です。`--specificity-output-root` のデフォルトは `../4_Specificity_Lat/outputs` です。
