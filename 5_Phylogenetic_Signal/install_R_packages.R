#!/usr/bin/env Rscript
packages <- c("ape", "phangorn", "phytools")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
message("Required R packages are available: ", paste(packages, collapse = ", "))
