#!/usr/bin/env Rscript

required <- c("bipartite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}
for (package in required) {
  cat(package, as.character(packageVersion(package)), "\n")
}
