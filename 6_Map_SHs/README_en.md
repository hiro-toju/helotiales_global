# 6_Map_SHs

This standalone tool draws one PDF occurrence map per abundant SH using the root unique-occurrence table produced by `1_GlobalFungi`.

It extracts only the old `*_selected_SH_spatial_maps/` functionality that was previously embedded in `2_Specificity`. It does not run randomizations, d′/2DP analyses, or FASTA extraction.

## Basic command

```bash
cd "[CURRENT DIRECTORY]"
MPLCONFIGDIR=/tmp/mplconfig python3 map_abundant_shs.py \
  --globalfungi-results ../1_GlobalFungi/results_helotiales \
  --output-dir outputs \
  --sh-selection rank \
  --sh-rank-max 50 \
  --figure-font-scale 1.5
```

This writes the top 50 SH maps to `outputs/sh50_maps/selected_SH_spatial_maps/`.

## Main options

- `--sh-selection rank`: select SHs by unique-occurrence rank.
- `--sh-rank-max 50`: number of top SHs to map.
- `--sh-selection min-count --sh-min-count 100`: map every SH with at least the requested number of unique occurrences.
- `--sh-selection list --selected-sh-file selected_shs.tsv`: map only explicitly listed SHs.
- `--keep-uncertain-hosts`: keep occurrences with ambiguous plant labels. Without this option, uncertain host labels are removed to mimic the old `2_Specificity` map output.
- `--figure-font-scale`: multiplier for figure text size.
- `--point-size`: map point size.

## Outputs

- `selected_SH_spatial_maps/*.pdf`: one geographic distribution map per selected SH
- `sh_unique_occurrence_ranking.tsv`: ranking table for all SHs
- `selected_shs.tsv`: SHs selected for mapping
- `excluded_occurrences.tsv`: occurrences excluded from mapping
- `run_summary.tsv`: input rows, excluded rows, mapped SHs, and runtime
- `command.txt`: executed command
- `package_versions.tsv`: Python and package versions
- `analysis.log`: compact workflow log

## Count unit

SH ranking is based on the number of unique occurrences defined by `sh_id + latitude + longitude + host_candidate`. This corresponds to the sampling-unit-like occurrence unit used in the current workflow.

