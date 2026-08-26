#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Usage: specificity_sample_unit_label_shuffle.R OCCURRENCE_TSV OUTPUT_PREFIX N SEED COLUMN_ROLE")
}

input_path <- args[[1L]]
output_prefix <- args[[2L]]
n_randomizations <- as.integer(args[[3L]])
seed <- as.integer(args[[4L]])
column_role <- args[[5L]]
if (!column_role %in% c("plant", "continent")) {
  stop("COLUMN_ROLE must be 'plant' or 'continent'")
}

if (!requireNamespace("bipartite", quietly = TRUE)) {
  stop("R package 'bipartite' is required")
}
if (is.na(n_randomizations) || n_randomizations < 2L) stop("N must be >= 2")

occ <- read.delim(
  input_path, check.names = FALSE, quote = "", comment.char = "",
  stringsAsFactors = FALSE, fileEncoding = "UTF-8"
)
required <- c(
  "occurrence_id", "sampling_unit_id", "spatial_block", "row_label",
  "col_label", "row_order", "col_order", "sample_ids"
)
if (!all(required %in% names(occ))) {
  stop("Occurrence input requires: ", paste(required, collapse = ", "))
}
if (nrow(occ) < 1L || anyNA(occ[, required])) stop("Occurrence input is empty or contains missing labels")
if (anyDuplicated(occ$occurrence_id)) stop("occurrence_id must be unique")
if (anyDuplicated(occ[, c("sampling_unit_id", "row_label")])) {
  stop("Each fungus must occur at most once in a sampling unit")
}

row_labels <- unique(occ[order(occ$row_order), c("row_label", "row_order")])$row_label
col_labels <- unique(occ[order(occ$col_order), c("col_label", "col_order")])$col_label
if (length(row_labels) < 2L || length(col_labels) < 2L) {
  stop("Randomization requires at least two fungal and two column labels")
}

matrix_from_labels <- function(column_labels) {
  unclass(table(
    factor(occ$row_label, levels = row_labels),
    factor(column_labels, levels = col_labels)
  ))
}

safe_dfun <- function(x) {
  row_ok <- rowSums(x) > 0
  col_ok <- colSums(x) > 0
  answer <- rep(NA_real_, nrow(x))
  if (sum(row_ok) == 0L || sum(col_ok) < 2L) return(answer)
  answer[row_ok] <- as.numeric(suppressWarnings(
    bipartite::dfun(x[row_ok, col_ok, drop = FALSE])$dprime
  ))
  answer
}

new_stats <- function(shape) {
  list(n = array(0L, dim = shape), mean = array(0, dim = shape), m2 = array(0, dim = shape))
}
update_stats <- function(stats, values) {
  valid <- is.finite(values)
  stats$n[valid] <- stats$n[valid] + 1L
  delta <- values[valid] - stats$mean[valid]
  stats$mean[valid] <- stats$mean[valid] + delta / stats$n[valid]
  stats$m2[valid] <- stats$m2[valid] + delta * (values[valid] - stats$mean[valid])
  stats
}
finish_stats <- function(stats) {
  sd <- array(NA_real_, dim = dim(stats$mean))
  enough <- stats$n > 1L
  sd[enough] <- sqrt(stats$m2[enough] / (stats$n[enough] - 1L))
  mean <- stats$mean
  mean[stats$n == 0L] <- NA_real_
  list(mean = mean, sd = sd, n = stats$n)
}
z_score <- function(observed, mean, sd) {
  z <- rep(NA_real_, length(observed)); dim(z) <- dim(observed)
  ok <- is.finite(observed) & is.finite(mean) & is.finite(sd) & sd > 0
  z[ok] <- (observed[ok] - mean[ok]) / sd[ok]
  z
}
empirical_tail <- function(exceed, valid_n) {
  p <- rep(NA_real_, length(exceed)); dim(p) <- dim(exceed)
  ok <- valid_n > 0L
  p[ok] <- (exceed[ok] + 1) / (valid_n[ok] + 1)
  p
}
empirical_two_sided <- function(ge, le, valid_n) {
  upper <- empirical_tail(ge, valid_n)
  lower <- empirical_tail(le, valid_n)
  p <- pmin(1, 2 * pmin(as.vector(upper), as.vector(lower)))
  dim(p) <- dim(valid_n)
  p
}
bh <- function(p) {
  out <- rep(NA_real_, length(p)); valid <- is.finite(p)
  out[valid] <- p.adjust(p[valid], method = "BH")
  out
}

split_sample_ids <- function(x) {
  if (length(x) == 0L || is.na(x) || !nzchar(trimws(x))) return(character(0))
  values <- unlist(strsplit(x, "[,;|]+"), use.names = FALSE)
  values <- trimws(values)
  values <- values[nzchar(values)]
  values <- values[!tolower(values) %in% c("unknown", "unidentified", "unresolved", "na", "nan", "none")]
  unique(values)
}

combine_sample_ids <- function(x) {
  values <- unique(unlist(lapply(x, split_sample_ids), use.names = FALSE))
  values <- values[nzchar(values)]
  paste(sort(values), collapse = "|")
}

build_permutation_plan <- function(indices, sample_sets) {
  n <- length(indices)
  incompatible <- matrix(FALSE, nrow = n, ncol = n)
  sample_to_positions <- new.env(hash = TRUE, parent = emptyenv())
  for (position in seq_len(n)) {
    for (sample_id in sample_sets[[indices[[position]]]]) {
      current <- if (exists(sample_id, sample_to_positions, inherits = FALSE)) {
        get(sample_id, sample_to_positions, inherits = FALSE)
      } else {
        integer(0)
      }
      assign(sample_id, c(current, position), sample_to_positions)
    }
  }
  for (sample_id in ls(sample_to_positions)) {
    positions <- get(sample_id, sample_to_positions, inherits = FALSE)
    incompatible[positions, positions] <- TRUE
  }
  if (any(rowSums(!incompatible) == 0L)) {
    stop("At least one sampling unit has no sample_id-compatible source")
  }
  if (any(colSums(!incompatible) == 0L)) {
    stop("At least one sampling unit cannot be used as a sample_id-compatible source")
  }
  list(indices = indices, incompatible = incompatible)
}

greedy_compatible_sources <- function(plan, max_attempts = 100L) {
  n <- length(plan$indices)
  if (n < 2L) return(seq_len(n))
  compatible <- !plan$incompatible
  for (attempt in seq_len(max_attempts)) {
    destination_order <- sample(seq_len(n))
    destination_order <- destination_order[order(rowSums(compatible)[destination_order])]
    assignment <- rep(NA_integer_, n)
    used <- rep(FALSE, n)
    failed <- FALSE
    for (dest in destination_order) {
      candidates <- which(compatible[dest, ] & !used)
      if (!length(candidates)) {
        failed <- TRUE
        break
      }
      chosen <- sample(candidates, 1L)
      assignment[dest] <- chosen
      used[chosen] <- TRUE
    }
    if (!failed && all(!is.na(assignment))) return(assignment)
  }
  stop(
    "Could not construct a sample_id-constrained label permutation after ",
    max_attempts, " attempts"
  )
}

fast_compatible_sources <- function(plan, max_attempts = 50L, max_repair_factor = 20L) {
  n <- length(plan$indices)
  if (n < 2L) return(seq_len(n))
  incompatible <- plan$incompatible
  if (!any(incompatible)) return(sample.int(n))
  destinations <- seq_len(n)
  max_repair_swaps <- max(1000L, as.integer(max_repair_factor) * n)
  for (attempt in seq_len(max_attempts)) {
    assignment <- sample.int(n)
    bad <- which(incompatible[cbind(destinations, assignment)])
    if (!length(bad)) return(assignment)
    for (repair_step in seq_len(max_repair_swaps)) {
      dest <- sample(bad, 1L)
      source <- assignment[[dest]]
      partners <- which(!incompatible[dest, assignment] & !incompatible[, source])
      partners <- partners[partners != dest]
      if (!length(partners)) break
      partner <- sample(partners, 1L)
      assignment[c(dest, partner)] <- assignment[c(partner, dest)]
      bad <- which(incompatible[cbind(destinations, assignment)])
      if (!length(bad)) return(assignment)
    }
  }
  greedy_compatible_sources(plan)
}

unit_table <- aggregate(
  sample_ids ~ sampling_unit_id + spatial_block + col_label,
  data = occ,
  FUN = combine_sample_ids
)
if (anyDuplicated(unit_table$sampling_unit_id)) {
  stop("Each sampling_unit_id must have exactly one spatial block and column label")
}
unit_sample_sets <- lapply(unit_table$sample_ids, split_sample_ids)
global_unit_indices <- seq_len(nrow(unit_table))
global_permutation_plan <- build_permutation_plan(global_unit_indices, unit_sample_sets)
observed <- matrix_from_labels(occ$col_label)
storage.mode(observed) <- "integer"
if (sum(observed) != nrow(occ)) stop("Observed matrix does not contain every occurrence")

observed_row <- safe_dfun(observed)
observed_col <- safe_dfun(t(observed))
row_stats <- new_stats(length(observed_row))
col_stats <- new_stats(length(observed_col))
cell_stats <- new_stats(dim(observed))
row_ge <- integer(length(observed_row)); row_le <- integer(length(observed_row))
col_ge <- integer(length(observed_col)); col_le <- integer(length(observed_col))
cell_ge <- matrix(0L, nrow(observed), ncol(observed))
cell_le <- matrix(0L, nrow(observed), ncol(observed))

set.seed(seed)
progress_step <- max(1L, n_randomizations %/% 10L)
loop_started_at <- Sys.time()
for (iteration in seq_len(n_randomizations)) {
  shuffled <- occ$col_label
  shuffled_units <- unit_table$col_label
  if (length(shuffled_units) > 1L) {
    source_order <- fast_compatible_sources(global_permutation_plan)
    shuffled_units <- shuffled_units[source_order]
  }
  shuffled <- shuffled_units[match(occ$sampling_unit_id, unit_table$sampling_unit_id)]
  randomized <- matrix_from_labels(shuffled)
  storage.mode(randomized) <- "integer"
  if (sum(randomized) != nrow(occ)) stop("Occurrence total changed after shuffling")

  randomized_row <- safe_dfun(randomized)
  randomized_col <- safe_dfun(t(randomized))
  row_stats <- update_stats(row_stats, randomized_row)
  col_stats <- update_stats(col_stats, randomized_col)
  cell_stats <- update_stats(cell_stats, randomized)
  row_valid <- is.finite(randomized_row) & is.finite(observed_row)
  col_valid <- is.finite(randomized_col) & is.finite(observed_col)
  row_ge[row_valid] <- row_ge[row_valid] + (randomized_row[row_valid] >= observed_row[row_valid])
  row_le[row_valid] <- row_le[row_valid] + (randomized_row[row_valid] <= observed_row[row_valid])
  col_ge[col_valid] <- col_ge[col_valid] + (randomized_col[col_valid] >= observed_col[col_valid])
  col_le[col_valid] <- col_le[col_valid] + (randomized_col[col_valid] <= observed_col[col_valid])
  cell_ge <- cell_ge + (randomized >= observed)
  cell_le <- cell_le + (randomized <= observed)
  if (iteration %% progress_step == 0L || iteration == n_randomizations) {
    elapsed_seconds <- as.numeric(difftime(Sys.time(), loop_started_at, units = "secs"))
    eta_seconds <- if (iteration > 0L) elapsed_seconds * (n_randomizations - iteration) / iteration else NA_real_
    message(
      "R global sample-unit ", column_role, "-label shuffle: ",
      iteration, "/", n_randomizations,
      " elapsed=", round(elapsed_seconds, 1), "s",
      " eta=", round(eta_seconds, 1), "s"
    )
  }
}

row_done <- finish_stats(row_stats); col_done <- finish_stats(col_stats); cell_done <- finish_stats(cell_stats)
row_z <- z_score(observed_row, row_done$mean, row_done$sd)
col_z <- z_score(observed_col, col_done$mean, col_done$sd)
cell_z <- z_score(observed, cell_done$mean, cell_done$sd)
row_upper <- empirical_tail(row_ge, row_done$n); row_lower <- empirical_tail(row_le, row_done$n)
col_upper <- empirical_tail(col_ge, col_done$n); col_lower <- empirical_tail(col_le, col_done$n)
row_directional <- ifelse(row_z < 0, row_lower, row_upper)
col_directional <- ifelse(col_z < 0, col_lower, col_upper)
row_p <- pmin(1, 2 * row_directional)
col_p <- pmin(1, 2 * col_directional)
cell_p <- empirical_two_sided(cell_ge, cell_le, cell_done$n)

row_table <- data.frame(
  label = row_labels, dprime_observed = observed_row,
  dprime_null_mean = row_done$mean, dprime_null_sd = row_done$sd,
  z_standardized_dprime = as.vector(row_z), p_randomization_upper = row_upper,
  p_randomization_lower = row_lower, p_randomization_directional_tail = row_directional,
  p_randomization_two_sided = row_p, fdr_bh_randomization_two_sided = bh(row_p),
  n_valid_randomizations = row_done$n, stringsAsFactors = FALSE
)
col_table <- data.frame(
  label = col_labels, dprime_observed = observed_col,
  dprime_null_mean = col_done$mean, dprime_null_sd = col_done$sd,
  z_standardized_dprime = as.vector(col_z), p_randomization_upper = col_upper,
  p_randomization_lower = col_lower, p_randomization_directional_tail = col_directional,
  p_randomization_two_sided = col_p, fdr_bh_randomization_two_sided = bh(col_p),
  n_valid_randomizations = col_done$n, stringsAsFactors = FALSE
)
indices <- expand.grid(row_index = seq_len(nrow(observed)), col_index = seq_len(ncol(observed)))
cell_p_vector <- cell_p[cbind(indices$row_index, indices$col_index)]
cell_table <- data.frame(
  row_label = row_labels[indices$row_index], col_label = col_labels[indices$col_index],
  observed_count = observed[cbind(indices$row_index, indices$col_index)],
  null_mean = cell_done$mean[cbind(indices$row_index, indices$col_index)],
  null_sd = cell_done$sd[cbind(indices$row_index, indices$col_index)],
  n_valid_randomizations = cell_done$n[cbind(indices$row_index, indices$col_index)],
  z_standardized_2dp = cell_z[cbind(indices$row_index, indices$col_index)],
  p_randomization_two_sided = cell_p_vector,
  fdr_bh_randomization_two_sided = bh(cell_p_vector), stringsAsFactors = FALSE
)

units_with_sample_ids <- vapply(unit_sample_sets, length, integer(1)) > 0L

constraints <- data.frame(
  model = "global_sampling_unit_label_shuffle",
  randomization_unit = "latitude_longitude_host_sampling_unit",
  shuffled_field = paste0(column_role, "_label"),
  specificity_domain = column_role,
  spatial_block = "none",
  n_binary_fungus_occurrences = nrow(occ), n_sampling_units = nrow(unit_table),
  n_sampling_units_with_sample_ids = sum(units_with_sample_ids),
  n_continents_in_input = length(unique(unit_table$spatial_block)),
  n_blocks = 1L,
  sampling_without_replacement = TRUE,
  forbids_label_exchange_between_units_sharing_sample_id = TRUE,
  n_forbidden_source_destination_pairs = sum(global_permutation_plan$incompatible),
  constrained_permutation_algorithm = "precomputed_incompatibility_fast_repair_with_greedy_fallback",
  constrained_permutation_fast_repair_max_attempts = 50L,
  constrained_permutation_greedy_fallback_max_attempts = 100L,
  preserves_fungal_cooccurrence_within_sampling_unit = TRUE,
  preserves_global_sampling_unit_total = TRUE,
  preserves_global_fungal_label_counts = TRUE,
  preserves_global_column_label_counts = TRUE,
  rebuilds_interaction_matrix_each_iteration = TRUE,
  n_randomizations = n_randomizations, seed = seed,
  stringsAsFactors = FALSE
)

write_table <- function(x, suffix) {
  write.table(
    x, paste0(output_prefix, suffix), sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE, na = "NA", fileEncoding = "UTF-8"
  )
}
write_table(row_table, "__row.tsv")
write_table(col_table, "__col.tsv")
write_table(cell_table, "__cell.tsv")
write_table(constraints, "__constraints.tsv")
