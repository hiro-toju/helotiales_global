# 4_Specificity_Lat: latitude-blocked randomization

This step applies the sampling-unit label shuffle with signed latitude bands instead of continent blocks. Northern and southern latitudes are kept separate; the default band width is 20 degrees.

Plant labels are not exchanged between sampling units sharing any `sample_id` / `sample_ids` value.

The constrained shuffle uses a faster implementation: sample-id-incompatible pairs are precomputed once, then only invalid assignments are repaired after an ordinary shuffle.

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

Fungal-genus-level example:

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

Pass this run name to `5_Phylogenetic_Signal --specificity-run-name`.
