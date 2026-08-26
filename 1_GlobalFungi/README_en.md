# 1_GlobalFungi: create Helotiales root occurrences from GlobalFungi

In version 2, root occurrences passed to downstream analyses are restricted to medium-or-higher host assignments. Output filenames do not include confidence-specific strings.

Main outputs:

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

`--jobs 0` uses all available CPU cores.
