# 7_Data_Overview

This auxiliary workflow creates overview figures for the Helotiales occurrence data.

It outputs world maps of sampling locations and continent-composition PDF plots for four nested data definitions:

1. All GlobalFungi occurrences.
2. All occurrences containing Helotiales.
3. All Helotiales occurrences from root samples, i.e. `sample_type == root`.
4. Root unique occurrences with reliable host-plant information.

For the first dataset, all GlobalFungi occurrences are counted as positive cells in the sample-by-SH abundance matrix in `GlobalFungi_5_SH_abundance_ITS1_ITS2.txt.gz`. The map shows sample coordinates for samples in which at least one SH was detected.

By default, the fourth dataset uses the cleaned table from `0_Data_Property`, in which 198 unique occurrences flagged during previous preprocessing are removed. If that file is not available, the script falls back to `1_GlobalFungi/results_helotiales/helotiales_root_occurrences_with_hosts.tsv`.

## Basic command

```bash
cd "[CURRENT DIRECTORY]"
python3 data_overview.py \
  --globalfungi-data ../1_GlobalFungi/data \
  --globalfungi-results ../1_GlobalFungi/results_helotiales \
  --output-dir outputs \
  --figure-font-scale 1.5
```

## Main outputs

- `outputs/data_overview_summary.tsv`  
  Numbers of records, records with coordinates, sampling locations, and continents for each data definition.

- `outputs/continent_breakdown.tsv`  
  Continent-level counts and proportions for each data definition.

- `outputs/figures/*_world_map.pdf`  
  World maps of sampling locations.

- `outputs/figures/*_continent_proportions.pdf`  
  Continent-composition plots. Continent order and colors are fixed across all figures, and continent names and percentages are shown directly inside or around the plots.

- `outputs/continent_color_key.tsv`  
  Fixed continent display order and colors.

- `outputs/run_configuration.tsv`, `runtime_versions.tsv`, `run_timing.tsv`, `workflow.log`  
  Reproducibility records for command settings, runtime versions, timing, and logs.
