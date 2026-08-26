#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
cd "$SCRIPT_DIR"
python3 phylogeny_workflow.py   --data-dir data   --specificity-run-name sh150_family30_10000_lat20_minblock100   --specificity-fasta-subset its1_and_its2   --target-taxon Helotiales   --outgroup-required-shared-rank class   --minimum-outgroup-genus-sh-count 3   --outgroup-minimum-median-identity 75   --outgroup-maximum-median-identity 85   --outgroup-target-median-identity 80   --outgroup-maximum-identity-to-any-target 90   --outgroup-minimum-median-query-coverage 75   --outgroup-minimum-query-hit-fraction 0.80   --trait-column host_z_standardized_dprime_sampling_unit_label_shuffle   --minimum-length 400   --outgroup-minimum-length 460   --terminal-trim-minimum-occupancy 1.0   --trimal-gap-threshold 0.50   --trimal-minimum-retained-percent 0   --bootstrap 500   --signal-randomizations 10000   --jobs 0   --iqtree-thread-mode auto   --figure-font-scale 1.2
