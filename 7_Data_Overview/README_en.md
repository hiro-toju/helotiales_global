# 7_Data_Overview

This auxiliary workflow creates overview figures for the Helotiales occurrence data.

It outputs world maps of sampling locations and continent-composition PDF plots for three nested data definitions:

1. All occurrences containing Helotiales.
2. All Helotiales occurrences from root samples, i.e. `sample_type == root`.
3. Root unique occurrences with reliable host-plant information.

By default, the third dataset uses the cleaned table from `0_Data_Property`, in which 198 unique occurrences flagged during previous preprocessing are removed. If that file is not available, the script falls back to `1_GlobalFungi/results_helotiales/helotiales_root_occurrences_with_hosts.tsv`.

## Basic command

```bash
cd "[CURRENT DIRECTORY]"
python3 data_overview.py \
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
  Continent-composition plots.

- `outputs/run_configuration.tsv`, `runtime_versions.tsv`, `run_timing.tsv`, `workflow.log`  
  Reproducibility records for command settings, runtime versions, timing, and logs.

