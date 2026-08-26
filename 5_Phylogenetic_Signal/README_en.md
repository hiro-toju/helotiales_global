# 5_Phylogenetic_Signal: phylogenetic signal analysis

Run this step last. By default, it automatically reads the integrated host-specificity summary, 2DP table, and FASTA from the matching `4_Specificity_Lat` run folder.

To use `3_Specificity_SpatialBlocks` results, explicitly set `--specificity-output-root ../3_Specificity_SpatialBlocks/outputs`.

## Example using 4_Specificity_Lat results

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

## Example using 3_Specificity_SpatialBlocks results

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

`--minimum-length` filters target sequences, whereas `--outgroup-minimum-length` filters UNITE outgroup candidates. The default `--specificity-output-root` is `../4_Specificity_Lat/outputs`.
