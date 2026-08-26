# 0_Data_Property

This step summarizes the basic properties of the cleaned Helotiales unique-occurrence data after applying the exclusions reported by `4_Specificity_Lat`.

## Inputs

The usual input is a completed `4_Specificity_Lat` run.

- `4_Specificity_Lat/outputs/<run_name>/data/helotiales_root_occurrences_with_hosts.tsv`
- `4_Specificity_Lat/outputs/<run_name>/data/plant_taxonomy.tsv`
- `4_Specificity_Lat/outputs/<run_name>/helotiales_root_occurrences_with_hosts/helotiales_root_occurrences_with_hosts_excluded_unique_occurrences.tsv`
- `4_Specificity_Lat/outputs/<run_name>/helotiales_root_occurrences_with_hosts/helotiales_root_occurrences_with_hosts_excluded_plant_columns.tsv`

## Basic command

Example for an SH-level run:

```bash
cd "[CURRENT DIRECTORY]"
MPLCONFIGDIR=/tmp/mplconfig python3 data_property_summary.py \
  --specificity-run-name sh150_family30_10000_lat20_minblock100 \
  --plot-top-n 15 \
  --top-sh-bar-n 50 \
  --figure-font-scale 1.2
```

Example for a fungal-genus-level run:

```bash
cd "[CURRENT DIRECTORY]"
MPLCONFIGDIR=/tmp/mplconfig python3 data_property_summary.py \
  --specificity-run-name genus150_family30_10000_lat20_minblock100 \
  --plot-top-n 15 \
  --top-sh-bar-n 50 \
  --figure-font-scale 1.2
```

`data_property_summary.py` switches inputs by run name, so the same Python script is used for both SH-level and fungal-genus-level runs.

## Main options

| Option | Description |
|---|---|
| `--specificity-run-name` | Name of the `4_Specificity_Lat` run to summarize |
| `--input` | Explicit occurrence TSV |
| `--specificity-result-dir` | Explicit folder containing exclusion reports |
| `--plant-taxonomy` | Explicit plant taxonomy TSV |
| `--output-dir` | Output directory. Defaults to `outputs/<run_name>` |
| `--plot-top-n` | Number of top categories shown individually in summary plots |
| `--top-sh-bar-n` | Number of top fungal SHs shown in the frequency bar plot. Default: 50 |
| `--figure-font-scale` | Font-size multiplier for figures |

## Outputs

Outputs are written under `outputs/<run_name>/`.

- Cleaned occurrence table after exclusions
- Count summaries for fungal SH, fungal genus, fungal family, plant genus, plant family, plant order, country, and continent
- PDF plots showing unique-occurrence composition for each category
- PDF bar plot for the most frequent fungal SHs
- Run command, timing, and Python/package version files

## Note

This is a descriptive data-property step. It does not replace the core downstream specificity inputs; it documents what remains after preprocessing filters.
