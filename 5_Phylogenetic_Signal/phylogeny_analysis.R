#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 11L || ((length(args) - 11L) %% 5L) != 0L) {
  stop("Usage: phylogeny_analysis.R ALIGNMENT ML_TREE TRAITS TWO_DP_TABLE OUTPUT OUTGROUP BOOTSTRAP SIGNAL_RANDOMIZATIONS SEED NJ_MODEL FONT_SCALE [EXTRA_TRAITS EXTRA_2DP TRAIT_SET TRAIT_TITLE LEGEND_TITLE ...]")
}
alignment_path <- args[[1]]; ml_path <- args[[2]]; trait_path <- args[[3]]
two_dp_path <- args[[4]]; output_dir <- args[[5]]; outgroup <- args[[6]]
bootstrap_n <- as.integer(args[[7]]); signal_n <- as.integer(args[[8]])
seed <- as.integer(args[[9]]); nj_model <- args[[10]]; font_scale <- as.numeric(args[[11]])
extra_args <- if (length(args) > 11L) args[12:length(args)] else character()

required <- c("ape", "phangorn", "phytools")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R package(s): ", paste(missing, collapse = ", "))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)

alignment <- ape::read.dna(alignment_path, format = "fasta", as.character = FALSE)
label_traits_df <- read.delim(trait_path, check.names = FALSE, stringsAsFactors = FALSE)

read_trait_set <- function(trait_file, two_dp_file, trait_set, trait_title, legend_title) {
  traits_df <- read.delim(trait_file, check.names = FALSE, stringsAsFactors = FALSE)
  two_dp <- read.delim(two_dp_file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!("sh_id" %in% names(two_dp))) stop("2DP table lacks sh_id: ", two_dp_file)
  two_dp <- two_dp[!duplicated(two_dp$sh_id), , drop = FALSE]
  list(
    trait_file = trait_file, two_dp_file = two_dp_file, trait_set = trait_set,
    trait_title = trait_title, legend_title = legend_title,
    traits_df = traits_df,
    traits = setNames(as.numeric(traits_df$trait_value), traits_df$tip_label),
    two_dp = two_dp
  )
}

trait_sets <- list(read_trait_set(
  trait_path, two_dp_path, "host", "host-specificity z-d'", "z-standardized d'"
))
if (length(extra_args)) {
  for (start in seq(1L, length(extra_args), by = 5L)) {
    trait_sets[[length(trait_sets) + 1L]] <- read_trait_set(
      extra_args[[start]], extra_args[[start + 1L]], extra_args[[start + 2L]],
      extra_args[[start + 3L]], extra_args[[start + 4L]]
    )
  }
}

specificity_z_palette <- grDevices::colorRampPalette(c(
  "#053061", "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
  "#F7F7F7",
  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"
))(201)

specificity_heatmap_z_limit <- function(values) {
  finite <- abs(as.numeric(values[is.finite(values)]))
  if (!length(finite)) return(3.0)
  min(max(3.0, as.numeric(stats::quantile(finite, 0.98, names = FALSE, type = 7))), 10.0)
}

make_symmetric_z_scale <- function(values, metric) {
  limit <- specificity_heatmap_z_limit(values)
  list(
    metric = metric,
    palette = specificity_z_palette,
    limits = c(-limit, limit),
    display_limit = limit
  )
}

dprime_scale <- make_symmetric_z_scale(
  unlist(lapply(trait_sets, function(set) set$traits), use.names = FALSE),
  "z_standardized_dprime"
)
two_dp_scale <- make_symmetric_z_scale(
  unlist(lapply(trait_sets, function(set) {
    two_dp <- set$two_dp[, setdiff(names(set$two_dp), "sh_id"), drop = FALSE]
    as.numeric(as.matrix(two_dp))
  }), use.names = FALSE),
  "z_standardized_2dp"
)
write.table(
  data.frame(
    metric = c(dprime_scale$metric, two_dp_scale$metric),
    color_map = "RdBu_r",
    low_color = "#053061",
    midpoint_value = 0,
    midpoint_color = "#F7F7F7",
    high_color = "#67001F",
    display_min = c(dprime_scale$limits[[1]], two_dp_scale$limits[[1]]),
    display_max = c(dprime_scale$limits[[2]], two_dp_scale$limits[[2]]),
    scale_rule = "global across trait-tree PDFs; symmetric around zero; limit = min(max(3, 98th percentile of |z|), 10)",
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "trait_tree_color_scales.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

root_tree <- function(tree) {
  if (!(outgroup %in% tree$tip.label)) stop("Outgroup is absent from tree: ", outgroup)
  original_distances <- ape::cophenetic.phylo(tree)
  outgroup_tip <- match(outgroup, tree$tip.label)
  original_edge_row <- which(tree$edge[, 2] == outgroup_tip)
  if (length(original_edge_row) != 1L) stop("Could not identify the outgroup pendant edge")
  original_pendant_length <- tree$edge.length[original_edge_row]
  rooted <- ape::root(tree, outgroup = outgroup, resolve.root = TRUE, edgelabel = TRUE)
  root_node <- ape::Ntip(rooted) + 1L
  root_rows <- which(rooted$edge[, 1] == root_node)
  outgroup_row <- which(rooted$edge[, 2] == match(outgroup, rooted$tip.label))
  other_row <- setdiff(root_rows, outgroup_row)
  if (length(root_rows) != 2L || length(outgroup_row) != 1L || length(other_row) != 1L) {
    stop("Outgroup rooting did not produce the expected outgroup-versus-ingroup bifurcation")
  }
  # Place the root at the midpoint of the original outgroup pendant edge. Splitting
  # that edge equally preserves every pairwise patristic distance.
  rooted$edge.length[c(outgroup_row, other_row)] <- original_pendant_length / 2
  rooted_distances <- ape::cophenetic.phylo(rooted)
  common <- intersect(rownames(original_distances), rownames(rooted_distances))
  distance_error <- max(abs(original_distances[common, common] - rooted_distances[common, common]))
  attr(rooted, "rooting_diagnostics") <- data.frame(
    outgroup = outgroup,
    original_outgroup_pendant_length = original_pendant_length,
    rooted_outgroup_half_edge_length = rooted$edge.length[outgroup_row],
    rooted_ingroup_half_edge_length = rooted$edge.length[other_row],
    maximum_pairwise_distance_error = distance_error,
    root_degree = length(root_rows),
    n_outgroup_tips = 1L,
    n_ingroup_tips = ape::Ntip(rooted) - 1L,
    whole_ingroup_bootstrap = NA_character_,
    whole_ingroup_bootstrap_note = "Not estimable with one outgroup: the one-tip-versus-rest split is trivial in every unrooted bootstrap tree.",
    stringsAsFactors = FALSE
  )
  rooted
}


clear_trivial_root_support <- function(tree) {
  if (is.null(tree$node.label)) return(tree)
  root_node <- ape::Ntip(tree) + 1L
  root_children <- tree$edge[tree$edge[, 1] == root_node, 2]
  outgroup_tip <- match(outgroup, tree$tip.label)
  ingroup_child <- setdiff(root_children, outgroup_tip)
  nodes <- c(root_node, ingroup_child[ingroup_child > ape::Ntip(tree)])
  tree$node.label[nodes - ape::Ntip(tree)] <- NA_character_
  tree
}

ml_unrooted <- ape::read.tree(ml_path)
alignment_labels <- rownames(alignment)
if (!setequal(ml_unrooted$tip.label, alignment_labels)) {
  stop("ML tree tips do not match the trimmed alignment; remove/replace stale IQ-TREE outputs and rerun")
}
ml_full <- root_tree(ml_unrooted)
ml_rooting_diagnostics <- attr(ml_full, "rooting_diagnostics")
ml_full <- clear_trivial_root_support(ml_full)
ape::write.tree(ml_full, file.path(output_dir, "07_maximum_likelihood_rooted.treefile"))

nj_fun <- function(x) {
  distance <- ape::dist.dna(x, model = nj_model, pairwise.deletion = TRUE)
  ape::nj(distance)
}
nj_full <- root_tree(nj_fun(alignment))
nj_rooting_diagnostics <- attr(nj_full, "rooting_diagnostics")
nj_boot <- ape::boot.phylo(
  nj_full, alignment, FUN = function(x) root_tree(nj_fun(x)), B = bootstrap_n, quiet = TRUE
)
nj_full$node.label <- sprintf("%.0f", 100 * nj_boot / bootstrap_n)
nj_full <- clear_trivial_root_support(nj_full)
ape::write.tree(nj_full, file.path(output_dir, "08_neighbor_joining_rooted.treefile"))
rooting_diagnostics <- rbind(
  cbind(tree_method = "maximum_likelihood", ml_rooting_diagnostics),
  cbind(tree_method = "neighbor_joining", nj_rooting_diagnostics)
)
write.table(rooting_diagnostics, file.path(output_dir, "outgroup_rooting_diagnostics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(
  data.frame(node = seq_along(nj_boot) + ape::Ntip(nj_full), bootstrap_percent = 100 * nj_boot / bootstrap_n),
  file.path(output_dir, "neighbor_joining_bootstrap.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)

prepare_trait <- function(full_tree, trait_values) {
  ingroup <- ape::drop.tip(full_tree, outgroup)
  finite <- is.finite(trait_values)
  trait_values <- trait_values[finite]
  common <- intersect(ingroup$tip.label, names(trait_values))
  ingroup <- ape::drop.tip(ingroup, setdiff(ingroup$tip.label, common))
  x <- trait_values[ingroup$tip.label]
  if (length(x) < 4L || any(!is.finite(x))) stop("At least four finite ingroup trait values are required")
  if (!ape::is.binary(ingroup)) ingroup <- ape::multi2di(ingroup, random = FALSE)
  ingroup$edge.length[!is.finite(ingroup$edge.length) | ingroup$edge.length <= 0] <- 1e-8
  list(tree = ingroup, trait = x)
}

signal_one <- function(full_tree, method_name, trait_values, trait_name) {
  prepared <- prepare_trait(full_tree, trait_values); tree <- prepared$tree; x <- prepared$trait
  k <- phytools::phylosig(tree, x, method = "K", test = TRUE, nsim = signal_n)
  lambda <- phytools::phylosig(tree, x, method = "lambda", test = TRUE)
  data.frame(
    tree_method = method_name, n_tips = length(x), trait = trait_name,
    blomberg_K = unname(k$K), blomberg_K_P = unname(k$P), blomberg_K_randomizations = signal_n,
    pagel_lambda = unname(lambda$lambda), pagel_lambda_logL = unname(lambda$logL),
    pagel_lambda_logL0 = unname(lambda$logL0), pagel_lambda_P = unname(lambda$P),
    stringsAsFactors = FALSE
  )
}

safe_trait_label <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

trait_label_column_name <- function(trait_set) {
  if (identical(trait_set, "host")) "host_family"
  else if (identical(trait_set, "continent")) "continent"
  else paste0(safe_trait_label(trait_set), "_label")
}

trait_set_fdr_stem <- function(trait_set) {
  if (identical(trait_set, "host")) "phylogenetic_conservatism_host_family"
  else paste0("phylogenetic_conservatism_", safe_trait_label(trait_set))
}

write_phylogenetic_signal_fdr_outputs <- function(signals, trait_set, labels_in_order) {
  if (!length(labels_in_order)) return(invisible(NULL))

  label_col <- trait_label_column_name(trait_set)
  safe_labels <- vapply(labels_in_order, safe_trait_label, character(1))
  trait_names <- paste0(trait_set, "_z_standardized_2dp_", safe_labels)
  label_lookup <- setNames(labels_in_order, trait_names)
  order_lookup <- setNames(seq_along(trait_names), trait_names)

  fdr_df <- signals[signals$trait %in% trait_names, , drop = FALSE]
  if (!nrow(fdr_df)) return(invisible(NULL))
  fdr_df[[label_col]] <- unname(label_lookup[fdr_df$trait])
  fdr_df$trait_order <- unname(order_lookup[fdr_df$trait])
  fdr_df$blomberg_K_fdr_bh <- NA_real_
  fdr_df$pagel_lambda_fdr_bh <- NA_real_

  for (method in unique(fdr_df$tree_method)) {
    idx <- which(fdr_df$tree_method == method)
    fdr_df$blomberg_K_fdr_bh[idx] <- p.adjust(fdr_df$blomberg_K_P[idx], method = "BH")
    fdr_df$pagel_lambda_fdr_bh[idx] <- p.adjust(fdr_df$pagel_lambda_P[idx], method = "BH")
  }
  fdr_df$fdr_method <- "Benjamini-Hochberg"
  fdr_df$fdr_scope <- paste0(
    "2DP traits within each tree_method for ", trait_set,
    "; Blomberg_K_P and Pagel_lambda_P corrected separately"
  )
  fdr_df <- fdr_df[order(fdr_df$tree_method, fdr_df$trait_order), , drop = FALSE]

  stem <- trait_set_fdr_stem(trait_set)
  combined_path <- file.path(output_dir, paste0(stem, "_fdr_bh.tsv"))
  write.table(fdr_df, combined_path, sep = "\t", quote = FALSE, row.names = FALSE)

  blomberg_cols <- c(
    "tree_method", label_col, "trait_order", "n_tips", "trait",
    "blomberg_K", "blomberg_K_P", "blomberg_K_fdr_bh",
    "blomberg_K_randomizations", "fdr_method", "fdr_scope"
  )
  pagel_cols <- c(
    "tree_method", label_col, "trait_order", "n_tips", "trait",
    "pagel_lambda", "pagel_lambda_logL", "pagel_lambda_logL0",
    "pagel_lambda_P", "pagel_lambda_fdr_bh", "fdr_method", "fdr_scope"
  )
  write.table(
    fdr_df[, intersect(blomberg_cols, names(fdr_df)), drop = FALSE],
    file.path(output_dir, paste0(stem, "_blomberg_K_fdr_bh.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  write.table(
    fdr_df[, intersect(pagel_cols, names(fdr_df)), drop = FALSE],
    file.path(output_dir, paste0(stem, "_pagel_lambda_fdr_bh.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  metrics <- c(
    "blomberg_K", "blomberg_K_P", "blomberg_K_fdr_bh",
    "pagel_lambda", "pagel_lambda_P", "pagel_lambda_fdr_bh"
  )
  wide_rows <- list()
  row_i <- 1L
  for (method in unique(fdr_df$tree_method)) {
    method_df <- fdr_df[fdr_df$tree_method == method, , drop = FALSE]
    method_df <- method_df[match(labels_in_order, method_df[[label_col]]), , drop = FALSE]
    for (metric in metrics) {
      values <- as.list(method_df[[metric]])
      names(values) <- labels_in_order
      wide_rows[[row_i]] <- c(
        list(tree_method = method, metric = metric, fdr_method = if (grepl("_fdr_bh$", metric)) "Benjamini-Hochberg" else "none"),
        values
      )
      row_i <- row_i + 1L
    }
  }
  wide <- as.data.frame(
    do.call(rbind, lapply(wide_rows, function(x) unlist(x, use.names = FALSE))),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(wide) <- names(wide_rows[[1]])
  write.table(
    wide,
    file.path(output_dir, paste0(stem, "_fdr_bh_by_", label_col, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  invisible(fdr_df)
}

plotmath_quote <- function(x) gsub('"', '\\\\"', x, fixed = TRUE)

join_plotmath_parts <- function(parts) paste(parts[nzchar(parts)], collapse = "~")

is_higher_taxon_label <- function(x) {
  grepl(
    "(aceae|ales|mycetes|mycotina|mycota|fungi|incertae|unidentified|unknown|^NA$)$",
    x, ignore.case = TRUE
  )
}

is_sp_epithet <- function(x) grepl("^sp\\.?$", x, ignore.case = TRUE)

tip_label_expression <- function(tip, show_count = TRUE) {
  if (tip == outgroup) {
    body <- sub("^OUTGROUP_", "", tip)
    hit <- regexec("^(SH[0-9.]+[A-Za-z]*)_(.+)_([^_]+)$", body)
    parts <- regmatches(body, hit)[[1]]
    if (length(parts) == 4L) {
      first_name <- parts[[3]]
      second_name <- parts[[4]]
      if (is_sp_epithet(second_name)) second_name <- "sp."
      first_part <- if (is_higher_taxon_label(first_name)) {
        sprintf('plain("%s")', plotmath_quote(first_name))
      } else {
        sprintf('italic("%s")', plotmath_quote(first_name))
      }
      second_part <- if (is_sp_epithet(second_name)) {
        'plain("sp.")'
      } else {
        sprintf('italic("%s")', plotmath_quote(second_name))
      }
      return(paste(
        sprintf('bold("OUTGROUP | %s")', plotmath_quote(parts[[2]])),
        first_part,
        second_part, sep = "~"
      ))
    }
    return(sprintf('bold("OUTGROUP | %s")', plotmath_quote(body)))
  }
  row <- label_traits_df[match(tip, label_traits_df$tip_label), , drop = FALSE]
  sid <- if (nrow(row) && "sh_id" %in% names(row)) row$sh_id[[1]] else tip
  count <- if (nrow(row) && "n_unique_occurrences" %in% names(row)) row$n_unique_occurrences[[1]] else {
    hit <- regmatches(if (nrow(row)) row$display_label[[1]] else tip,
                      regexpr("(?<=n=)[0-9]+", if (nrow(row)) row$display_label[[1]] else tip, perl = TRUE))
    ifelse(length(hit) && nzchar(hit), hit, "0")
  }
  species <- if (nrow(row) && "species" %in% names(row)) row$species[[1]] else {
    label <- if (nrow(row)) row$display_label[[1]] else tip
    sub("_\\[n=.*$", "", sub(paste0("^", sid, "_"), "", label))
  }
  species <- gsub("_", " ", species)
  genus <- if (nrow(row) && "genus" %in% names(row)) row$genus[[1]] else strsplit(species, " +")[[1]][[1]]
  words <- strsplit(trimws(species), " +")[[1]]
  if (length(words) >= 2L && is_sp_epithet(words[[2]])) words[[2]] <- "sp."
  epithet <- if (length(words) >= 2L) words[[2]] else ""
  valid_genus <- nzchar(genus) && !is_higher_taxon_label(genus) && !is_higher_taxon_label(words[[1]])
  sid_part <- sprintf('plain("%s")', plotmath_quote(sid))
  suffix_part <- if (show_count) sprintf('plain("[n=%s]")', plotmath_quote(as.character(count))) else ""
  if (!valid_genus) {
    plain_name <- trimws(paste(words, collapse = " "))
    return(join_plotmath_parts(c(sid_part, sprintf('plain("%s")', plotmath_quote(plain_name)), suffix_part)))
  }
  genus_part <- sprintf('italic("%s")', plotmath_quote(genus))
  if (!nzchar(epithet) || is_sp_epithet(epithet)) {
    return(join_plotmath_parts(c(sid_part, genus_part, 'plain("sp.")', suffix_part)))
  }
  epithet_part <- sprintf('italic("%s")', plotmath_quote(epithet))
  join_plotmath_parts(c(sid_part, genus_part, epithet_part, suffix_part))
}

plot_trait_tree <- function(full_tree, trait_values, trait_title, legend_title, color_scale,
                            method_name, stem, support_mode = c("all", "gt70")) {
  support_mode <- match.arg(support_mode)
  # The first constrained tip is plotted at the bottom, placing the outgroup and
  # ancestral root below the ingroup rather than above it.
  full_tree <- ape::rotateConstr(full_tree, c(outgroup, setdiff(full_tree$tip.label, outgroup)))
  prepared <- prepare_trait(full_tree, trait_values); ingroup <- prepared$tree; x <- prepared$trait
  ancestral <- phytools::fastAnc(ingroup, x)
  palette <- color_scale$palette
  limits <- color_scale$limits
  value_color <- function(value) {
    value <- as.numeric(value)
    colors <- rep("#C4C4C4", length(value))
    ok <- is.finite(value)
    if (any(ok)) {
      index <- 1L + round(200 * (value[ok] - limits[[1]]) / diff(limits))
      colors[ok] <- palette[pmax(1L, pmin(201L, index))]
    }
    colors
  }
  descendant_tips <- function(tree, node) {
    tips <- integer(); frontier <- node
    while (length(frontier)) {
      children <- tree$edge[tree$edge[, 1] %in% frontier, 2]
      tips <- c(tips, children[children <= ape::Ntip(tree)])
      frontier <- children[children > ape::Ntip(tree)]
    }
    tree$tip.label[sort(unique(tips))]
  }
  node_values <- rep(NA_real_, ape::Ntip(full_tree) + full_tree$Nnode)
  names(node_values) <- seq_along(node_values)
  for (tip in intersect(full_tree$tip.label, names(x))) {
    node_values[match(tip, full_tree$tip.label)] <- x[[tip]]
  }
  ingroup_root <- ape::Ntip(ingroup) + 1L
  root_value <- unname(ancestral[as.character(ingroup_root)])
  full_root <- ape::Ntip(full_tree) + 1L
  node_values[full_root] <- root_value
  for (node in seq.int(ape::Ntip(full_tree) + 1L, ape::Ntip(full_tree) + full_tree$Nnode)) {
    tips <- intersect(setdiff(descendant_tips(full_tree, node), outgroup), names(x))
    if (!length(tips)) next
    if (length(tips) == 1L) {
      node_values[node] <- x[[tips]]
    } else {
      mrca <- ape::getMRCA(ingroup, tips)
      node_values[node] <- if (is.null(mrca)) root_value else unname(ancestral[as.character(mrca)])
    }
  }
  edge_values <- vapply(seq_len(nrow(full_tree$edge)), function(i) {
    parent <- node_values[full_tree$edge[i, 1]]; child <- node_values[full_tree$edge[i, 2]]
    if (is.finite(parent) && is.finite(child)) mean(c(parent, child))
    else if (is.finite(child)) child else if (is.finite(parent)) parent else root_value
  }, numeric(1))
  edge_colors <- value_color(edge_values)
  outgroup_tip <- match(outgroup, full_tree$tip.label)
  edge_colors[full_tree$edge[, 2] == outgroup_tip] <- "#4D4D4D"
  plot_tree <- full_tree
  plot_tree$edge.length <- plot_tree$edge.length * 0.72
  original_labels <- plot_tree$tip.label
  show_tip_counts <- identical(support_mode, "all")
  label_expressions <- as.expression(lapply(
    vapply(original_labels, tip_label_expression, character(1), show_count = show_tip_counts),
    function(label) parse(text = label)[[1]]
  ))
  tip_colors <- ifelse(original_labels == outgroup, "#4D4D4D", "#111111")
  draw <- function(device) {
    old <- par(family = "Arial", mar = c(3.5, 0.6, 4.4, 1.1), xpd = NA); on.exit(par(old), add = TRUE)
    tree_depth <- max(ape::node.depth.edgelength(plot_tree))
    n_tips <- ape::Ntip(plot_tree)
    choose_tip_cex <- function() {
      lower <- 0.35 * font_scale
      upper <- 0.96 * font_scale
      available_height <- max(par("pin")[[2]] * 0.90, 1)
      for (i in seq_len(24)) {
        mid <- (lower + upper) / 2
        needed <- n_tips * strheight("Ag", units = "inches", cex = mid) * 1.08
        if (needed <= available_height) lower <- mid else upper <- mid
      }
      lower
    }
    tip_cex <- 0.81 * choose_tip_cex()
    label_inches <- max(strwidth(label_expressions, units = "inches", cex = tip_cex))
    label_space <- label_inches * tree_depth / max(par("pin")[[1]] - 0.5, 1)
    ape::plot.phylo(
      plot_tree, type = "phylogram", show.tip.label = FALSE,
      edge.color = edge_colors, edge.width = 3,
      x.lim = c(0, tree_depth + label_space + 0.03 * tree_depth),
      y.lim = c(-1.55, ape::Ntip(plot_tree) + 0.9), no.margin = FALSE
    )
    plot_env <- get(".PlotPhyloEnv", envir = asNamespace("ape"))
    last_plot <- get("last_plot.phylo", envir = plot_env)
    text(last_plot$xx[seq_len(ape::Ntip(plot_tree))] + 0.012 * tree_depth,
         last_plot$yy[seq_len(ape::Ntip(plot_tree))], labels = label_expressions,
         adj = c(0, 0.5), cex = tip_cex, col = tip_colors)
    if (!is.null(plot_tree$node.label)) {
      support <- suppressWarnings(as.numeric(plot_tree$node.label))
      root_label_index <- (ape::Ntip(plot_tree) + 1L) - ape::Ntip(plot_tree)
      support[root_label_index] <- NA_real_
      show_support <- is.finite(support) & (support_mode == "all" | support >= 70)
      labels <- ifelse(show_support, sprintf("%.0f", support), "")
      internal_nodes <- seq.int(ape::Ntip(plot_tree) + 1L, ape::Ntip(plot_tree) + plot_tree$Nnode)
      text(
        last_plot$xx[internal_nodes] - 0.006 * tree_depth,
        last_plot$yy[internal_nodes] + 0.14,
        labels = labels,
        adj = c(1, 0),
        cex = 0.405 * font_scale
      )
    }
    usr <- par("usr"); bar_x <- seq(usr[[1]] + 0.03 * diff(usr[1:2]), usr[[1]] + 0.28 * diff(usr[1:2]), length.out = 202)
    bar_y0 <- -1.15; bar_y1 <- -0.72
    rect(bar_x[-202], bar_y0, bar_x[-1], bar_y1, col = palette, border = NA)
    text(bar_x[[1]], bar_y0 - 0.18, labels = format(limits[[1]], digits = 3), adj = c(0, 1), cex = 0.55 * font_scale)
    text(mean(range(bar_x)), bar_y0 - 0.18, labels = "0", adj = c(0.5, 1), cex = 0.55 * font_scale)
    text(bar_x[[202]], bar_y0 - 0.18, labels = format(limits[[2]], digits = 3), adj = c(1, 1), cex = 0.55 * font_scale)
    text(mean(range(bar_x)), bar_y1 + 0.12, labels = legend_title, cex = 0.58 * font_scale)
    mtext(paste0(method_name, ": ", trait_title), side = 3, line = 2.1, cex = 1.0 * font_scale, font = 2)
    mtext("Outgroup-rooted at the midpoint of its pendant edge; outgroup excluded from trait reconstruction",
          side = 3, line = 0.8, cex = 0.68 * font_scale)
  }
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = 10.5, height = 11.5, family = "Arial", bg = "white")
  draw("pdf"); dev.off()
}

write_trait_set_outputs <- function(set) {
  traits_df <- set$traits_df
  traits <- set$traits
  two_dp <- set$two_dp
  trait_set <- set$trait_set
  safe_set <- gsub("[^A-Za-z0-9._-]+", "_", trait_set)
  tree_dir <- file.path(output_dir, "trait_trees")
  support_folder <- function(support_mode) {
    if (support_mode == "all") "bootstrap_all" else "bootstrap_70_or_higher"
  }
  dprime_dir_for <- function(support_mode) {
    base <- file.path(tree_dir, support_folder(support_mode))
    if (trait_set == "host") file.path(base, "z_standardized_dprime")
    else file.path(base, paste0(safe_set, "_z_standardized_dprime"))
  }
  two_dp_dir_for <- function(support_mode) {
    base <- file.path(tree_dir, support_folder(support_mode))
    if (trait_set == "host") file.path(base, "z_standardized_2dp")
    else file.path(base, paste0(safe_set, "_z_standardized_2dp"))
  }
  conservatism_file <- if (trait_set == "host") {
    file.path(output_dir, "phylogenetic_conservatism.tsv")
  } else {
    file.path(output_dir, paste0("phylogenetic_conservatism_", safe_set, ".tsv"))
  }

  for (support_mode in c("all", "gt70")) {
    dprime_dir <- dprime_dir_for(support_mode)
    plot_trait_tree(ml_full, traits, set$trait_title, set$legend_title, dprime_scale,
                    "Maximum likelihood", file.path(dprime_dir, "maximum_likelihood"), support_mode)
    plot_trait_tree(nj_full, traits, set$trait_title, set$legend_title, dprime_scale,
                    "Neighbor joining", file.path(dprime_dir, "neighbor_joining"), support_mode)
  }

  tip_by_sh <- setNames(traits_df$tip_label, traits_df$sh_id)
  plant_columns <- setdiff(names(two_dp), "sh_id")
  two_dp_signals <- list()
  for (plant in plant_columns) {
    values <- suppressWarnings(as.numeric(two_dp[[plant]]))
    names(values) <- unname(tip_by_sh[as.character(two_dp$sh_id)])
    values <- values[!is.na(names(values)) & is.finite(values)]
    safe_plant <- gsub("[^A-Za-z0-9._-]+", "_", plant)
    for (support_mode in c("all", "gt70")) {
      two_dp_dir <- two_dp_dir_for(support_mode)
      plot_trait_tree(ml_full, values, paste0(plant, " z-standardized 2DP"), "z-standardized 2DP", two_dp_scale,
                      "Maximum likelihood", file.path(two_dp_dir, paste0(safe_plant, "_maximum_likelihood")), support_mode)
      plot_trait_tree(nj_full, values, paste0(plant, " z-standardized 2DP"), "z-standardized 2DP", two_dp_scale,
                      "Neighbor joining", file.path(two_dp_dir, paste0(safe_plant, "_neighbor_joining")), support_mode)
    }
    two_dp_signals[[safe_plant]] <- rbind(
      signal_one(ml_full, "maximum_likelihood", values, paste0(trait_set, "_z_standardized_2dp_", safe_plant)),
      signal_one(nj_full, "neighbor_joining", values, paste0(trait_set, "_z_standardized_2dp_", safe_plant))
    )
  }
  dprime_signals <- rbind(
    signal_one(ml_full, "maximum_likelihood", traits, paste0(trait_set, "_z_standardized_dprime")),
    signal_one(nj_full, "neighbor_joining", traits, paste0(trait_set, "_z_standardized_dprime"))
  )
  signals <- if (length(two_dp_signals)) {
    do.call(rbind, c(list(dprime_signals), two_dp_signals))
  } else {
    dprime_signals
  }
  write.table(signals, conservatism_file, sep = "\t", quote = FALSE, row.names = FALSE)
  write_phylogenetic_signal_fdr_outputs(signals, trait_set, plant_columns)
  signals
}

all_signals <- do.call(rbind, lapply(trait_sets, write_trait_set_outputs))
write.table(all_signals, file.path(output_dir, "phylogenetic_conservatism_all_trait_sets.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

session <- capture.output(sessionInfo())
writeLines(session, file.path(output_dir, "R_sessionInfo.txt"))
