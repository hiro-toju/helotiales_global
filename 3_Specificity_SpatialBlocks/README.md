# 3_Specificity_SpatialBlocks

Continent-blocked host-specificity randomization. The workflow reads the medium-or-higher GlobalFungi output and uses neutral run-level filenames.

Plant labels are shuffled within continent blocks. Label exchange is forbidden between sampling units sharing any `sample_id` / `sample_ids` value.

The constrained shuffle uses a faster implementation: sample-id-incompatible pairs are precomputed once, then only invalid assignments are repaired after an ordinary shuffle.

```bash
cd "[CURRENT DIRECTORY]"
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

Fungal-genus-level example:

```bash
cd "[CURRENT DIRECTORY]"
MPLCONFIGDIR=/tmp/mplconfig python3 spatial_blocks_genus_workflow.py \
  --globalfungi-results ../1_GlobalFungi/results_helotiales \
  --output-dir outputs \
  --plant-rank family \
  --fungi-selection rank --fungi-rank-max 150 \
  --plant-selection rank --plant-rank-max 30 \
  --n-randomizations 10000 \
  --jobs 0 \
  --offline
```
