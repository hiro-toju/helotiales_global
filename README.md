# Analysis code for “Global distributions and host plant associations of the ecologically diverse ascomycete order Helotiales”

This repository contains reproducible Python/R analysis workflows associated with the manuscript “Global distributions and plant associations of the ecologically diverse ascomycete order Helotiales”.

The workflows analyze global occurrence records of root-associated Helotiales fungi using GlobalFungi and UNITE database files. They identify Helotiales Species Hypotheses (SHs), curate root-associated occurrence records with host-plant information, define unique occurrences and sampling units, summarize fungal and plant taxonomic composition, quantify host-plant and geographic specificity using sampling-unit label-shuffle randomizations, perform spatial-block and latitude-block randomization analyses, generate distribution maps for abundant SHs, and test phylogenetic signal in host-association specificity using ITS/rRNA sequence alignments and phylogenetic trees.

## Repository structure

The repository is organized as a stepwise workflow:

- `0_Data_Property/`  
  Summarizes the curated Helotiales occurrence dataset and visualizes the composition of fungal taxa, plant taxa, countries, and continents.

- `1_GlobalFungi/`  
  Processes GlobalFungi and UNITE input files, extracts Helotiales root-associated records, assigns host-plant information, defines unique occurrences, and prepares downstream input files.

- `2_Specificity/`  
  Performs host-plant and continent specificity analyses without spatial blocking.

- `3_Specificity_SpatialBlocks/`  
  Performs host-plant specificity analyses using continent-based spatial blocks.

- `4_Specificity_Lat/`  
  Performs host-plant specificity analyses using latitude-based spatial blocks.

- `5_Phylogenetic_Signal/`  
  Extracts SH sequences, selects an outgroup from UNITE, performs sequence alignment and phylogenetic inference, and tests phylogenetic signal in host-association specificity.

- `6_Map_SHs/`  
  Generates occurrence maps for abundant SHs.

## Required external data

Large primary database files from GlobalFungi and UNITE are not redistributed in this repository. To reproduce the analyses, download the following files from their original sources and place them in the indicated folders.

From GlobalFungi release 5.0 (https://globalfungi.com/), place the following files in `1_GlobalFungi/data/`:

- `GlobalFungi_5_SH_abundance_ITS1_ITS2.txt.gz`
- `GlobalFungi_5_sample_metadata.txt.gz`

From the UNITE database (https://unite.ut.ee/repository.php), place the following files in `1_GlobalFungi/data/`:

- `sh_general_release_s_04.04.2024.tgz`
- `sh_general_release_dynamic_04.04.2024.SHs.tax.bz2`

For the phylogenetic-signal workflow, the UNITE sequence archive is also required in `5_Phylogenetic_Signal/data/`:

- `sh_general_release_s_04.04.2024.tgz`

## Software environment

The workflows were developed and tested primarily on macOS. Installation notes and example commands assume a macOS environment with Homebrew-managed command-line tools unless otherwise stated. The scripts may also run on Linux-like environments if the required Python/R packages and external tools are installed and available on the system `PATH`, but non-macOS platforms were not the primary target of optimization.

The analyses use Python and R scripts together with external bioinformatics tools, including MAFFT, trimAl, BLAST+, and IQ-TREE/ModelFinder. Required R packages include `ape`, `phangorn`, `phytools`, and `bipartite`.

Each workflow records command-line settings, runtime information, software/package versions, and timing logs to support reproducibility.

## Notes on AI-assisted code development

Parts of the analysis code and documentation were drafted and refactored with assistance from ChatGPT/OpenAI Codex under the direction of the authors. The authors reviewed, tested, and revised the code and outputs, and take responsibility for the final implementation and scientific interpretation.
