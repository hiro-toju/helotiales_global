# 1_GlobalFungi：GlobalFungiからHelotiales root occurrenceを作成

第2版では、下流解析に渡すroot occurrenceを medium-or-higher host assignment のみに統一し、出力ファイル名にはconfidence表記を使いません。

主な出力：

```text
results_helotiales/helotiales_root_occurrences_with_hosts.tsv
results_helotiales/helotiales_root_host_sequences_all.fasta
results_helotiales/helotiales_root_host_sequences_its1_only.fasta
results_helotiales/helotiales_root_host_sequences_its2_only.fasta
results_helotiales/helotiales_root_host_sequences_its1_and_its2.fasta
```

```bash
cd "[CURRENT DIRECTORY]"
MPLCONFIGDIR=/tmp/mplconfig python3 helotiales_globalfungi_workflow.py \
  --occurrences data/GlobalFungi_5_SH_abundance_ITS1_ITS2.txt.gz \
  --metadata data/GlobalFungi_5_sample_metadata.txt.gz \
  --taxonomy data/sh_general_release_dynamic_04.04.2024.SHs.tax.bz2 \
  --fasta data/sh_general_release_s_04.04.2024.tgz \
  --target-taxon Helotiales \
  --figure-font-scale 1.5 \
  --output results_helotiales \
  --jobs 0
```

`--jobs 0` は利用可能CPU数を使います。
