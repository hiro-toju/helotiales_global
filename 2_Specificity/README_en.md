# 2_Specificity: unblocked randomization

This step reads the medium-or-higher output from `1_GlobalFungi/results_helotiales` and copies it into each run folder using the neutral name `helotiales_root_occurrences_with_hosts.tsv`.

Randomization shuffles plant labels and continent labels of sampling units. Label exchange is forbidden between sampling units sharing any `sample_id` / `sample_ids` value.

The constrained shuffle uses a faster implementation: sample-id-incompatible pairs are precomputed once, then only invalid assignments are repaired after an ordinary shuffle.

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

Example for a fungal-genus-level run:

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

`--jobs 0` uses all available CPU cores.
