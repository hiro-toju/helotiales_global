#!/usr/bin/env bash
set -euo pipefail

cd "[CURRENT DIRECTORY]"
python3 helotiales_globalfungi_workflow.py \
  --occurrences data/GlobalFungi_5_SH_abundance_ITS1_ITS2.txt.gz \
  --metadata data/GlobalFungi_5_sample_metadata.txt.gz \
  --taxonomy data/sh_general_release_dynamic_04.04.2024.SHs.tax.bz2 \
  --target-taxon Helotiales \
  --figure-font-scale 1.5 \
  --output results_helotiales
