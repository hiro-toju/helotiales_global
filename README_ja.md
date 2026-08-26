# Helotiales_Specificity_2

Helotiales root-associated fungal specificity workflow, version 2.

この第2版では、解析対象を medium-or-higher host assignment のみに統一します。第2版内部のrun用ファイル名では、旧版の `confidence_medium_or_high` のような文字列を使わず、`helotiales_root_occurrences_with_hosts.tsv` のような汎用名で扱います。

## ワークフロー順序

```text
0_Data_Property
1_GlobalFungi
2_Specificity
3_Specificity_SpatialBlocks
4_Specificity_Lat
5_Phylogenetic_Signal
```

系統解析は最後です。`5_Phylogenetic_Signal` は `4_Specificity_Lat` の出力runフォルダから、宿主特異性の統合表、2DP表、FASTAを自動的に探します。

`6_Map_SHs` は本流の後続解析ではなく、`1_GlobalFungi` のunique occurrence表からabundant SHの分布図だけを作る任意の補助ステップです。

## 第2版の主な仕様

- 下流解析では medium-or-higher の宿主判定のみを使います。
- `2_Specificity`, `3_Specificity_SpatialBlocks`, `4_Specificity_Lat` では、1_GlobalFungiのmedium-or-higher出力をコピーし、run内ではconfidence表記のないファイル名にリネームします。
- ランダマイゼーションでは、同じ `sample_id` または `sample_ids` を共有するsampling unit間で、植物ラベル・continentラベルを交換しません。
- IQ-TREE以外の並列化では、12を上限とする人工的な制約を入れていません。`--jobs 0` は利用可能CPU数を使います。
- 図はArial、PDF出力を基本とします。
- `outputs/` と `tmp/` は `.gitignore` で公開対象外です。

## 代表的な実行順

### 1. GlobalFungi

```bash
cd "[CURRENT DIRECTORY]/1_GlobalFungi"
MPLCONFIGDIR=/tmp/mplconfig python3 helotiales_globalfungi_workflow.py \
  --occurrences data/GlobalFungi_5_SH_abundance_ITS1_ITS2.txt.gz \
  --metadata data/GlobalFungi_5_sample_metadata.txt.gz \
  --taxonomy data/sh_general_release_dynamic_04.04.2024.SHs.tax.bz2 \
  --fasta data/sh_general_release_s_04.04.2024.tgz \
  --target-taxon Helotiales \
  --output results_helotiales \
  --jobs 0
```

### 2. blockなし specificity

```bash
cd "[CURRENT DIRECTORY]/2_Specificity"
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

### 3. continent spatial blocks

```bash
cd "[CURRENT DIRECTORY]/3_Specificity_SpatialBlocks"
MPLCONFIGDIR=/tmp/mplconfig python3 spatial_blocks_sh_workflow.py \
  --globalfungi-results ../1_GlobalFungi/results_helotiales \
  --output-dir outputs \
  --plant-rank family \
  --fungi-selection rank --fungi-rank-max 150 \
  --plant-selection rank --plant-rank-max 30 \
  --n-randomizations 10000 \
  --jobs 0 \
  --offline
```

### 4. latitude blocks

```bash
cd "[CURRENT DIRECTORY]/4_Specificity_Lat"
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

### 5. phylogenetic signal

`--specificity-run-name` には、4_Specificity_Latで作成されたrun名を指定します。

```bash
cd "[CURRENT DIRECTORY]/5_Phylogenetic_Signal"
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

### 任意. abundant SH maps

```bash
cd "[CURRENT DIRECTORY]/6_Map_SHs"
MPLCONFIGDIR=/tmp/mplconfig python3 map_abundant_shs.py \
  --globalfungi-results ../1_GlobalFungi/results_helotiales \
  --output-dir outputs \
  --sh-selection rank \
  --sh-rank-max 50 \
  --figure-font-scale 1.5
```
