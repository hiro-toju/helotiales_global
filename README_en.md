# Helotiales_Specificity_2

Helotiales root-associated fungal specificity workflow, version 2.

Version 2 analyzes only medium-or-higher host assignments. Inside v2 run folders, files are given neutral names such as `helotiales_root_occurrences_with_hosts.tsv`; old confidence-specific strings are not used in run-level filenames.

## Workflow order

```text
0_Data_Property
1_GlobalFungi
2_Specificity
3_Specificity_SpatialBlocks
4_Specificity_Lat
5_Phylogenetic_Signal
```

Phylogenetic analysis is the final step. `5_Phylogenetic_Signal` reads host-specificity summaries, 2DP tables, and FASTA files from the corresponding `4_Specificity_Lat` run folder.

`6_Map_SHs` is an optional side step, not a downstream analytical step. It creates abundant-SH distribution maps directly from the `1_GlobalFungi` unique-occurrence table.

## Main v2 rules

- Downstream analyses use only medium-or-higher host assignments.
- `2_Specificity`, `3_Specificity_SpatialBlocks`, and `4_Specificity_Lat` copy the medium-or-higher GlobalFungi outputs and rename them to neutral filenames inside each run folder.
- During randomization, labels are not exchanged between sampling units sharing any `sample_id` or `sample_ids` value.
- Outside IQ-TREE, `--jobs` has no artificial cap of 12. `--jobs 0` uses all available CPU cores.
- Figures use Arial and vector PDF output where possible.
- `outputs/` and `tmp/` are excluded by `.gitignore`.

## Optional abundant-SH maps

```bash
cd "[CURRENT DIRECTORY]/6_Map_SHs"
MPLCONFIGDIR=/tmp/mplconfig python3 map_abundant_shs.py \
  --globalfungi-results ../1_GlobalFungi/results_helotiales \
  --output-dir outputs \
  --sh-selection rank \
  --sh-rank-max 50 \
  --figure-font-scale 1.5
```

See `README_ja.md` for full example commands.
