#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(limma)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(openssl)
})

ROOT <- Sys.getenv("VS_INPUT_ROOT", unset = ".")
EXPORT_ROOT <- Sys.getenv("VS_OUTPUT_ROOT", unset = ".")
ANALYSIS <- file.path(EXPORT_ROOT, "outputs")
SCRIPTS <- file.path(EXPORT_ROOT, "scripts")
FIGURES <- file.path(EXPORT_ROOT, "figures")
SEED <- 20261015L
B <- 500L
VIF_LIMIT <- 10
CONDITION_LIMIT <- 1000
MIN_EFFECT_FOR_RATIO <- 0.1

dir.create(ANALYSIS, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES, recursive = TRUE, showWarnings = FALSE)

log_lines <- c(
  paste("bulk-model/effect-agreement independent refit started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste("R:", R.version.string),
  paste("Bootstrap seed:", SEED, "replicates:", B)
)
log_note <- function(...) {
  msg <- paste0(...)
  log_lines <<- c(log_lines, paste(format(Sys.time(), "%H:%M:%S"), msg))
  message(msg)
}

write_csv <- function(x, path) {
  data.table::fwrite(as.data.frame(x), path, na = "NA", quote = "auto")
}

read_expr <- function(path) {
  tab <- data.table::fread(path, check.names = FALSE, data.table = FALSE)
  genes <- as.character(tab[[1]])
  if (anyDuplicated(genes)) stop("Duplicate expression gene IDs: ", path)
  mat <- as.matrix(tab[, -1, drop = FALSE])
  storage.mode(mat) <- "double"
  rownames(mat) <- genes
  mat
}

read_csv <- function(path) {
  data.table::fread(path, check.names = FALSE, data.table = FALSE, encoding = "UTF-8")
}

sha256_file <- function(path) {
  con <- file(path, open = "rb")
  on.exit(close(con))
  as.character(openssl::sha256(con))
}

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop("Validation failed: ", message)
}

check_sample_alignment <- function(expr, meta, label) {
  if (!all(c("sample", "group") %in% names(meta))) stop(label, " metadata lacks sample/group")
  set_equal <- setequal(colnames(expr), meta$sample)
  meta_ordered <- meta[match(colnames(expr), meta$sample), , drop = FALSE]
  order_equal <- identical(as.character(meta_ordered$sample), colnames(expr))
  if (!set_equal || !order_equal || anyNA(meta_ordered$sample)) {
    stop(label, " expression/metadata sample alignment failed")
  }
  list(meta = meta_ordered, same_sample_set = set_equal, same_sample_order = order_equal)
}

fraction_read <- function(path, column_map, meta, label) {
  tab <- read_csv(path)
  if (!"sample" %in% names(tab)) stop(label, " fraction table lacks sample")
  if (!all(names(column_map) %in% names(tab))) {
    stop(label, " missing fraction columns: ", paste(setdiff(names(column_map), names(tab)), collapse = ", "))
  }
  if (!setequal(as.character(tab$sample), as.character(meta$sample))) {
    stop(label, " fraction sample set does not match discovery metadata")
  }
  tab <- tab[match(meta$sample, tab$sample), , drop = FALSE]
  frac <- as.matrix(tab[, names(column_map), drop = FALSE])
  storage.mode(frac) <- "double"
  colnames(frac) <- unname(column_map)
  rownames(frac) <- as.character(tab$sample)
  if (any(!is.finite(frac)) || any(frac < -1e-8)) stop(label, " fractions contain nonfinite/negative values")
  if (any(abs(rowSums(frac) - 1) > 1e-5)) stop(label, " fractions do not sum to one within 1e-5")
  list(table = tab, fraction = frac, sample_match = TRUE, row_sum_max_abs_error = max(abs(rowSums(frac) - 1)))
}

composition_pcs <- function(frac, n_pc = 5L) {
  frac <- as.matrix(frac)
  if (any(!is.finite(frac)) || any(frac < 0)) stop("Invalid fraction matrix in CLR-PCA")
  logged <- log(pmax(frac, 1e-6))
  clr <- logged - rowMeans(logged)
  pc <- stats::prcomp(clr, center = TRUE, scale. = FALSE)
  if (ncol(pc$x) < n_pc) stop("Requested more PCs than available")
  scores <- pc$x[, seq_len(n_pc), drop = FALSE]
  colnames(scores) <- paste0("PC", seq_len(n_pc))
  rownames(scores) <- rownames(frac)
  list(scores = scores, variance_fraction = pc$sdev^2 / sum(pc$sdev^2), clr = clr)
}

make_design <- function(meta, pcs = NULL) {
  dat <- data.frame(
    group = factor(as.character(meta$group), levels = c("Control", "Tumor")),
    batch = if ("batch" %in% names(meta)) factor(as.character(meta$batch)) else NULL
  )
  if (!is.null(pcs)) {
    pcs <- as.matrix(pcs)
    if (nrow(pcs) != nrow(dat)) stop("PC/sample count mismatch")
    for (j in seq_len(ncol(pcs))) dat[[colnames(pcs)[j]]] <- pcs[, j]
  }
  if ("batch" %in% names(meta)) {
    stats::model.matrix(stats::as.formula(paste("~ group + batch", if (!is.null(pcs)) paste0(" + ", paste(colnames(pcs), collapse = " + ")) else "")), data = dat)
  } else {
    stats::model.matrix(stats::as.formula(paste("~ group", if (!is.null(pcs)) paste0(" + ", paste(colnames(pcs), collapse = " + ")) else "")), data = dat)
  }
}

design_diagnostics <- function(design, meta, model_id, model_class, reference,
                               pc_n, analysis_set, baseline_vif = NA_real_,
                               original_main = FALSE) {
  rank <- qr(design)$rank
  full_rank <- rank == ncol(design)
  condition_approx <- if (full_rank) tryCatch(as.numeric(kappa(design)), error = function(e) Inf) else Inf
  condition_exact <- if (full_rank) tryCatch(as.numeric(kappa(design, exact = TRUE)), error = function(e) Inf) else Inf
  group_col <- grep("^groupTumor$", colnames(design))
  group_r2 <- NA_real_
  group_vif <- Inf
  if (length(group_col) == 1L) {
    g <- design[, group_col]
    other_cols <- setdiff(seq_len(ncol(design)), group_col)
    nuisance <- design[, other_cols, drop = FALSE]
    nuisance_fit <- tryCatch(stats::lm.fit(nuisance, g), error = function(e) NULL)
    if (!is.null(nuisance_fit) && sum((g - mean(g))^2) > 0 && all(is.finite(nuisance_fit$fitted.values))) {
      group_r2 <- 1 - sum(nuisance_fit$residuals^2) / sum((g - mean(g))^2)
      group_r2 <- min(1, max(0, group_r2))
      group_vif <- if (group_r2 < 1) 1 / (1 - group_r2) else Inf
    }
  }
  reasons <- character()
  if (!full_rank) reasons <- c(reasons, "rank deficient")
  if (!is.finite(condition_approx) || condition_approx > CONDITION_LIMIT) reasons <- c(reasons, "original approximate condition number > 1000/nonfinite")
  if (!is.finite(group_vif) || group_vif > VIF_LIMIT) reasons <- c(reasons, "group VIF > 10/nonfinite")
  data.frame(
    model_id = model_id,
    model_class = model_class,
    reference = reference,
    method = if (reference == "none") "none" else reference,
    pc_n = pc_n,
    pc_terms = if (pc_n == 0L) "none" else paste0("PC1-PC", pc_n),
    analysis_set = analysis_set,
    n = nrow(design),
    rank = rank,
    n_columns = ncol(design),
    residual_df = nrow(design) - rank,
    full_rank = full_rank,
    condition_number = condition_approx,
    original_approx_condition_number = condition_approx,
    exact_2norm_condition_number = condition_exact,
    condition_method = "gate: default approximate kappa(design), as in original script; separately report exact 2-norm kappa(design, exact=TRUE); unscaled design",
    group_R2_from_nuisance = group_r2,
    group_VIF = group_vif,
    baseline_group_VIF = baseline_vif,
    group_VIF_ratio_vs_baseline = if (is.finite(baseline_vif) && baseline_vif > 0) group_vif / baseline_vif else NA_real_,
    accepted_by_original_numeric_gate = full_rank && is.finite(condition_approx) && condition_approx <= CONDITION_LIMIT && is.finite(group_vif) && group_vif <= VIF_LIMIT,
    gate_failure_reason = if (!length(reasons)) "" else paste(reasons, collapse = "; "),
    original_main_MuSiC_PC1_5 = original_main,
    stringsAsFactors = FALSE
  )
}

fit_limma <- function(expr, meta, design, model_diag, model_id) {
  if (!isTRUE(model_diag$full_rank[1])) return(NULL)
  if (!"groupTumor" %in% colnames(design)) stop("Expected groupTumor coefficient in ", model_id)
  fit <- limma::lmFit(expr[, meta$sample, drop = FALSE], design)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  coef_name <- "groupTumor"
  coeff <- fit$coefficients[, coef_name]
  pval <- fit$p.value[, coef_name]
  fdr <- stats::p.adjust(pval, method = "BH")
  unscaled <- fit$stdev.unscaled[, coef_name]
  raw_se <- fit$sigma * unscaled
  moderated_se <- sqrt(fit$s2.post) * unscaled
  tab <- data.frame(
    gene = rownames(expr),
    model_id = model_id,
    logFC = as.numeric(coeff),
    P.Value = as.numeric(pval),
    adj.P.Val = as.numeric(fdr),
    residual_sigma = as.numeric(fit$sigma),
    stdev_unscaled = as.numeric(unscaled),
    raw_SE = as.numeric(raw_se),
    moderated_SE = as.numeric(moderated_se),
    moderated_t = as.numeric(fit$t[, coef_name]),
    posterior_variance = as.numeric(fit$s2.post),
    residual_df = as.numeric(fit$df.residual),
    total_df = as.numeric(fit$df.total),
    stringsAsFactors = FALSE
  )
  list(table = tab, fit = fit)
}

summarize_model <- function(model_diag, gene_table, baseline_table, baseline_sig) {
  base_sig_genes <- baseline_table$gene[baseline_sig]
  comp <- gene_table[match(base_sig_genes, gene_table$gene), , drop = FALSE]
  bs <- baseline_table[match(base_sig_genes, baseline_table$gene), , drop = FALSE]
  effect_sig <- !is.na(comp$adj.P.Val) & comp$adj.P.Val < 0.05 & abs(comp$logFC) > 1
  fdr_pass <- !is.na(comp$adj.P.Val) & comp$adj.P.Val < 0.05
  effect_pass <- abs(comp$logFC) > 1
  status <- ifelse(fdr_pass & effect_pass, "both thresholds retained",
                   ifelse(!fdr_pass & effect_pass, "only FDR threshold lost",
                          ifelse(fdr_pass & !effect_pass, "only effect threshold lost", "both thresholds lost")))
  base_status_counts <- table(factor(status, levels = c("both thresholds retained", "only FDR threshold lost", "only effect threshold lost", "both thresholds lost")))
  data.frame(
    model_id = model_diag$model_id,
    model_class = model_diag$model_class,
    reference = model_diag$reference,
    pc_n = model_diag$pc_n,
    analysis_set = model_diag$analysis_set,
    accepted_by_original_numeric_gate = model_diag$accepted_by_original_numeric_gate,
    original_main_MuSiC_PC1_5 = model_diag$original_main_MuSiC_PC1_5,
    baseline_fresh_DEG_n = length(base_sig_genes),
    adjusted_fdr05_logFC1_n_all_genes = sum(gene_table$adj.P.Val < 0.05 & abs(gene_table$logFC) > 1, na.rm = TRUE),
    both_thresholds_retained = unname(base_status_counts["both thresholds retained"]),
    only_FDR_threshold_lost = unname(base_status_counts["only FDR threshold lost"]),
    only_effect_threshold_lost = unname(base_status_counts["only effect threshold lost"]),
    both_thresholds_lost = unname(base_status_counts["both thresholds lost"]),
    sum_baseline_DEG_states = sum(base_status_counts),
    adjusted_only_threshold_calls_all_genes = sum(!gene_table$gene %in% base_sig_genes & gene_table$adj.P.Val < 0.05 & abs(gene_table$logFC) > 1, na.rm = TRUE),
    direction_reversals_among_fresh_baseline_DEGs = sum(sign(bs$logFC) != sign(comp$logFC), na.rm = TRUE),
    median_absolute_logFC_change_baseline_DEGs = stats::median(abs(comp$logFC - bs$logFC), na.rm = TRUE),
    median_raw_SE_ratio_baseline_DEGs = stats::median(comp$raw_SE / bs$raw_SE, na.rm = TRUE),
    median_moderated_SE_ratio_baseline_DEGs = stats::median(comp$moderated_SE / bs$moderated_SE, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

make_threshold_status <- function(gene_table, baseline_table, model_id, analysis_set) {
  base <- baseline_table
  base_sig <- !is.na(base$adj.P.Val) & base$adj.P.Val < 0.05 & abs(base$logFC) > 1
  base <- base[base_sig, c("gene", "logFC", "P.Value", "adj.P.Val", "raw_SE", "moderated_SE"), drop = FALSE]
  names(base)[-1] <- paste0("baseline_", names(base)[-1])
  comp <- gene_table[, c("gene", "logFC", "P.Value", "adj.P.Val", "raw_SE", "moderated_SE"), drop = FALSE]
  names(comp)[-1] <- paste0("scenario_", names(comp)[-1])
  out <- merge(base, comp, by = "gene", all.x = TRUE, sort = FALSE)
  fdr_pass <- !is.na(out$scenario_adj.P.Val) & out$scenario_adj.P.Val < 0.05
  effect_pass <- abs(out$scenario_logFC) > 1
  out$threshold_status <- ifelse(fdr_pass & effect_pass, "both thresholds retained",
                                 ifelse(!fdr_pass & effect_pass, "only FDR threshold lost",
                                        ifelse(fdr_pass & !effect_pass, "only effect threshold lost", "both thresholds lost")))
  out$direction_reversed <- sign(out$baseline_logFC) != sign(out$scenario_logFC)
  out$model_id <- model_id
  out$analysis_set <- analysis_set
  out$absolute_logFC_change <- abs(out$scenario_logFC - out$baseline_logFC)
  out$raw_SE_ratio <- out$scenario_raw_SE / out$baseline_raw_SE
  out$moderated_SE_ratio <- out$scenario_moderated_SE / out$baseline_moderated_SE
  out
}

range_overlap <- function(a, b, feature, group_name, reference) {
  a <- a[is.finite(a)]
  b <- b[is.finite(b)]
  if (!length(a) || !length(b)) return(data.frame(reference = reference, feature = feature, group1 = "Control", group2 = "Tumor", control_min = NA_real_, control_max = NA_real_, tumor_min = NA_real_, tumor_max = NA_real_, overlap_width = NA_real_, overlap_fraction_of_union_range = NA_real_))
  lo <- max(min(a), min(b))
  hi <- min(max(a), max(b))
  union_width <- max(max(a), max(b)) - min(min(a), min(b))
  data.frame(reference = reference, feature = feature, group1 = "Control", group2 = "Tumor",
             control_min = min(a), control_max = max(a), tumor_min = min(b), tumor_max = max(b),
             overlap_width = max(0, hi - lo),
             overlap_fraction_of_union_range = if (union_width > 0) max(0, hi - lo) / union_width else NA_real_)
}

distribution_rows <- function(values, meta, feature, reference) {
  do.call(rbind, lapply(c("Control", "Tumor"), function(g) {
    v <- values[meta$group == g]
    data.frame(reference = reference, feature = feature, group = g, n = sum(is.finite(v)),
               mean = mean(v, na.rm = TRUE), sd = stats::sd(v, na.rm = TRUE),
               median = stats::median(v, na.rm = TRUE), min = min(v, na.rm = TRUE), max = max(v, na.rm = TRUE),
               stringsAsFactors = FALSE)
  }))
}

effect_metrics <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- as.numeric(x[ok])
  y <- as.numeric(y[ok])
  n <- length(x)
  if (n < 3L) {
    return(c(n_pair = n, direction_n = 0, direction_match_n = NA_real_, direction_rate = NA_real_,
             Pearson = NA_real_, Spearman = NA_real_, CCC = NA_real_, mean_discovery = NA_real_,
             SD_discovery = NA_real_, mean_comparator = NA_real_, SD_comparator = NA_real_,
             mean_difference = NA_real_, OLS_slope = NA_real_, OLS_intercept = NA_real_))
  }
  nonzero <- x != 0 & y != 0
  direction_n <- sum(nonzero)
  direction_match_n <- if (direction_n) sum(sign(x[nonzero]) == sign(y[nonzero])) else NA_real_
  direction_rate <- if (direction_n) direction_match_n / direction_n else NA_real_
  varx_pop <- mean((x - mean(x))^2)
  vary_pop <- mean((y - mean(y))^2)
  cov_pop <- mean((x - mean(x)) * (y - mean(y)))
  pearson <- if (varx_pop > 0 && vary_pop > 0) stats::cor(x, y, method = "pearson") else NA_real_
  spearman <- if (length(unique(x)) > 1L && length(unique(y)) > 1L) suppressWarnings(stats::cor(x, y, method = "spearman")) else NA_real_
  ccc_denom <- varx_pop + vary_pop + (mean(x) - mean(y))^2
  ccc <- if (ccc_denom > 0) 2 * cov_pop / ccc_denom else NA_real_
  slope <- if (varx_pop > 0) cov_pop / varx_pop else NA_real_
  intercept <- if (is.finite(slope)) mean(y) - slope * mean(x) else NA_real_
  c(n_pair = n, direction_n = direction_n, direction_match_n = direction_match_n,
    direction_rate = direction_rate, Pearson = pearson, Spearman = spearman, CCC = ccc,
    mean_discovery = mean(x), SD_discovery = stats::sd(x),
    mean_comparator = mean(y), SD_comparator = stats::sd(y),
    mean_difference = mean(y - x), OLS_slope = slope, OLS_intercept = intercept)
}

fdr_rank_stat <- function(genes, fdr_x, fdr_y, comparison_id, scope) {
  ix <- match(genes, names(fdr_x))
  iy <- match(genes, names(fdr_y))
  x <- fdr_x[ix]
  y <- fdr_y[iy]
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  rank_cor <- if (length(x) >= 3L && length(unique(x)) > 1L && length(unique(y)) > 1L) {
    suppressWarnings(stats::cor(rank(x, ties.method = "average"), rank(y, ties.method = "average"), method = "spearman"))
  } else NA_real_
  data.frame(comparison_id = comparison_id, scope = scope,
             statistic = "Spearman correlation of full-fit BH-FDR ranks; significance ranking, not effect agreement",
             n_pair = length(x), FDR_rank_Spearman = rank_cor, stringsAsFactors = FALSE)
}

matrix_ols_group_coef <- function(expr, sample_ids, design, coefficient = "groupTumor") {
  if (qr(design)$rank != ncol(design)) return(NULL)
  if (!coefficient %in% colnames(design)) return(NULL)
  coef_matrix <- tryCatch(qr.solve(design, t(expr[, sample_ids, drop = FALSE])), error = function(e) NULL)
  if (is.null(coef_matrix)) return(NULL)
  stats::setNames(as.numeric(coef_matrix[coefficient, ]), rownames(expr))
}

external_group_fit <- function(expr, meta, contrast_id, tumor_samples, control_samples) {
  samples <- c(control_samples, tumor_samples)
  dat <- meta[match(samples, meta$sample), , drop = FALSE]
  if (anyNA(dat$sample)) stop("External sample metadata missing for ", contrast_id)
  dat$group <- factor(ifelse(dat$sample %in% control_samples, "Control", "Tumor"), levels = c("Control", "Tumor"))
  design <- stats::model.matrix(~ group, data = dat)
  if (qr(design)$rank != ncol(design)) stop("External group-only design rank deficient: ", contrast_id)
  fit <- limma::eBayes(limma::lmFit(expr[, samples, drop = FALSE], design), trend = TRUE, robust = TRUE)
  coeff <- fit$coefficients[, "groupTumor"]
  direct <- rowMeans(expr[, tumor_samples, drop = FALSE]) - rowMeans(expr[, control_samples, drop = FALSE])
  max_diff <- max(abs(coeff - direct), na.rm = TRUE)
  if (!is.finite(max_diff) || max_diff > 1e-10) stop("OLS group coefficient differs from direct group means: ", contrast_id, " diff=", max_diff)
  tab <- data.frame(
    gene = rownames(expr), contrast_id = contrast_id,
    n_control = length(control_samples), n_tumor = length(tumor_samples),
    logFC = as.numeric(coeff), P.Value = as.numeric(fit$p.value[, "groupTumor"]),
    adj.P.Val = stats::p.adjust(fit$p.value[, "groupTumor"], method = "BH"),
    residual_sigma = as.numeric(fit$sigma), stdev_unscaled = as.numeric(fit$stdev.unscaled[, "groupTumor"]),
    raw_SE = as.numeric(fit$sigma * fit$stdev.unscaled[, "groupTumor"]),
    moderated_SE = as.numeric(sqrt(fit$s2.post) * fit$stdev.unscaled[, "groupTumor"]),
    moderated_t = as.numeric(fit$t[, "groupTumor"]), posterior_variance = as.numeric(fit$s2.post),
    residual_df = as.numeric(fit$df.residual), total_df = as.numeric(fit$df.total),
    stringsAsFactors = FALSE
  )
  list(table = tab, max_ols_mean_difference = max_diff, n = nrow(dat), design_rank = qr(design)$rank, design_columns = ncol(design))
}

sample_stratified_draw <- function(meta, strata) {
  ix <- integer()
  draw <- list()
  for (nm in names(strata)) {
    stratum_index <- which(strata[[nm]])
    if (!length(stratum_index)) stop("Empty bootstrap stratum: ", nm)
    selected <- stratum_index[sample.int(length(stratum_index), length(stratum_index), replace = TRUE)]
    draw[[nm]] <- selected
    ix <- c(ix, selected)
  }
  list(index = ix, by_stratum = draw)
}

format_ids <- function(indices, sample_ids) paste(sample_ids[indices], collapse = ";")

emit_boot_metric_rows <- function(replicate, comparison_id, metrics, observed_gate_failure,
                                  fit_estimable, failure_reason = NULL) {
  names_keep <- c("direction_n", "direction_match_n", "direction_rate", "Pearson", "Spearman", "CCC", "mean_discovery",
                  "SD_discovery", "mean_comparator", "SD_comparator", "mean_difference",
                  "OLS_slope", "OLS_intercept")
  do.call(rbind, lapply(names_keep, function(nm) {
    val <- if (nm %in% names(metrics)) as.numeric(metrics[[nm]]) else NA_real_
    valid <- length(val) == 1L && is.finite(val)
    reason <- if (!is.null(failure_reason)) failure_reason else if (valid) "" else "metric not estimable in this draw (e.g. constant vector or no nonzero pairs)"
    data.frame(replicate = replicate, comparison_id = comparison_id, metric = nm,
               value = if (valid) val else NA_real_, valid = valid,
               failure_reason = reason, n_pair = if ("n_pair" %in% names(metrics)) as.integer(metrics[["n_pair"]]) else 0L,
               observed_composition_gate_failure = observed_gate_failure,
               fit_estimable = fit_estimable, stringsAsFactors = FALSE)
  }))
}

summarize_bootstrap <- function(draw_table, B, seed) {
  ids <- unique(draw_table$comparison_id)
  metrics <- unique(draw_table$metric)
  do.call(rbind, lapply(ids, function(cid) do.call(rbind, lapply(metrics, function(met) {
    rows <- draw_table[draw_table$comparison_id == cid & draw_table$metric == met, , drop = FALSE]
    v <- rows$value[rows$valid & is.finite(rows$value)]
    data.frame(comparison_id = cid, metric = met, bootstrap_seed = seed,
               requested_B = B, n_valid = length(v), n_invalid = B - length(v),
               valid_draw_fraction = length(v) / B,
               interval_scope = "Conditional on fixed reference fractions/annotations and fixed submitted 346-gene set; sample resampling only",
               percentile_2.5 = if (length(v)) unname(stats::quantile(v, 0.025, names = FALSE, type = 7)) else NA_real_,
               percentile_97.5 = if (length(v)) unname(stats::quantile(v, 0.975, names = FALSE, type = 7)) else NA_real_,
               interval_method = "percentile; finite pre-determined draws only; no redraw",
               stringsAsFactors = FALSE)
  }))))
}

write_plot_formats <- function(plot, stem, width_in, height_in) {
  if (!requireNamespace("svglite", quietly = TRUE) || !requireNamespace("ragg", quietly = TRUE)) {
    stop("Figure output requires svglite and ragg in the selected R runtime")
  }
  svg_path <- paste0(stem, ".svg")
  pdf_path <- paste0(stem, ".pdf")
  png_path <- paste0(stem, ".png")
  svglite::svglite(svg_path, width = width_in, height = height_in)
  print(plot)
  grDevices::dev.off()
  grDevices::cairo_pdf(pdf_path, width = width_in, height = height_in, family = "Arial")
  print(plot)
  grDevices::dev.off()
  ragg::agg_png(png_path, width = width_in, height = height_in, units = "in", res = 600)
  print(plot)
  grDevices::dev.off()
  invisible(c(svg_path, pdf_path, png_path))
}

theme_publication <- function() {
  ggplot2::theme_classic(base_size = 8, base_family = "Arial") +
    ggplot2::theme(
      axis.text = ggplot2::element_text(family = "Arial", size = 7),
      axis.title = ggplot2::element_text(family = "Arial", size = 8),
      strip.text = ggplot2::element_text(family = "Arial", size = 7, face = "bold"),
      plot.title = ggplot2::element_text(family = "Arial", size = 8, face = "bold"),
      plot.subtitle = ggplot2::element_text(family = "Arial", size = 6.5),
      legend.title = ggplot2::element_text(family = "Arial", size = 7),
      legend.text = ggplot2::element_text(family = "Arial", size = 6.5),
      panel.grid = ggplot2::element_blank()
    )
}

run_main <- function() {
  paths <- list(
    discovery_expr = file.path(ROOT, "data/processed/GSE39645_expr_gene.csv"),
    discovery_meta = file.path(ROOT, "data/processed/GSE39645_meta.csv"),
    external_expr = file.path(ROOT, "data/processed/GSE108524_expr_gene.csv"),
    external_meta = file.path(ROOT, "data/processed/GSE108524_meta.csv"),
    hubs = file.path(ANALYSIS, "locked_submitted_hub_set_346.csv"),
    source_samples = file.path(ANALYSIS, "control_sample_sources_9.csv"),
    subtype_evidence = file.path(ROOT, "results/nf2_sporadic/sample_label_evidence.csv"),
    music = file.path(ROOT, "results/deconv/cell_fractions_GSE39645_music.csv"),
    nnls = file.path(ROOT, "results/deconv/cell_fractions_GSE39645_nnls.csv"),
    gse216783 = file.path(ROOT, "results/gse216783/cell_fractions_GSE39645_gse216783ref.csv")
  )
  for (p in unlist(paths)) assert_true(file.exists(p), paste("input missing:", p))
  input_hash <- data.frame(path = unlist(paths), sha256 = vapply(unlist(paths), sha256_file, character(1)), stringsAsFactors = FALSE)
  write_csv(input_hash, file.path(ANALYSIS, "bulk_agreement_INPUT_SHA256.csv"))

  log_note("Reading processed matrices, locked hubs, metadata, and fixed reference fractions")
  expr <- read_expr(paths$discovery_expr)
  meta_raw <- read_csv(paths$discovery_meta)
  ext_expr <- read_expr(paths$external_expr)
  ext_meta <- read_csv(paths$external_meta)
  hubtab <- read_csv(paths$hubs)
  hubs <- as.character(hubtab$gene)
  assert_true(length(hubs) == 346L && !anyDuplicated(hubs), "locked identity has exactly 346 unique genes")
  assert_true(all(hubs %in% rownames(expr)) && all(hubs %in% rownames(ext_expr)), "locked genes measurable in both cohorts")
  assert_true(all(c("Control", "Tumor") %in% unique(meta_raw$group)), "discovery groups are Control/Tumor")
  aligned <- check_sample_alignment(expr, meta_raw, "GSE39645")
  meta <- aligned$meta
  meta$group <- factor(as.character(meta$group), levels = c("Control", "Tumor"))
  meta$batch <- factor(as.character(meta$batch))
  assert_true(nrow(meta) == 40L && sum(meta$group == "Tumor") == 31L && sum(meta$group == "Control") == 9L,
              "GSE39645 expected 31 Tumor and 9 controls")
  ext_aligned <- check_sample_alignment(ext_expr, ext_meta, "GSE108524")
  ext_meta <- ext_aligned$meta
  ext_meta$group <- factor(as.character(ext_meta$group), levels = c("Control", "Tumor"))
  assert_true(nrow(ext_meta) == 31L && sum(ext_meta$group == "Tumor") == 27L && sum(ext_meta$group == "Control") == 4L,
              "GSE108524 expected 27 Tumor and 4 controls")

  colmaps <- list(
    MuSiC = c("Schwann/肿瘤细胞" = "Schwann", "巨噬/髓系" = "Myeloid", "T细胞" = "T_cell", "成纤维" = "Fibroblast", "内皮" = "Endothelial", "周细胞" = "Pericyte"),
    NNLS = c("Schwann/肿瘤细胞" = "Schwann", "巨噬/髓系" = "Myeloid", "T细胞" = "T_cell", "成纤维" = "Fibroblast", "内皮" = "Endothelial", "周细胞" = "Pericyte"),
    GSE216783_reference = c("Schwann/tumour" = "Schwann", "Myeloid" = "Myeloid", "T cell" = "T_cell", "Fibroblast" = "Fibroblast", "Endothelial" = "Endothelial", "Pericyte" = "Pericyte")
  )
  fraction_files <- list(MuSiC = paths$music, NNLS = paths$nnls, GSE216783_reference = paths$gse216783)
  fractions <- lapply(names(fraction_files), function(nm) fraction_read(fraction_files[[nm]], colmaps[[nm]], meta, nm))
  names(fractions) <- names(fraction_files)

  source_tab <- read_csv(paths$source_samples)
  tissue_samples <- unique(c(as.character(meta$sample[meta$group == "Tumor"]),
    as.character(source_tab$sample[source_tab$source != "Human Schwann Cells"])))
  meta_tissue <- meta[meta$sample %in% tissue_samples, , drop = FALSE]
  expr_tissue <- expr[, meta_tissue$sample, drop = FALSE]
  assert_true(nrow(meta_tissue) == 39L && sum(meta_tissue$group == "Tumor") == 31L && sum(meta_tissue$group == "Control") == 8L,
              "tissue-only analysis is 31 Tumor vs 8 tissue controls")
  meta_tissue$group <- factor(as.character(meta_tissue$group), levels = c("Control", "Tumor"))
  meta_tissue$batch <- factor(as.character(meta_tissue$batch), levels = levels(meta$batch))

  # Fixed bulk-model model universe: one baseline, twelve whole-cohort reference/PC fits,
  # plus tissue-only baseline and MuSiC PC1-5 diagnostics.
  model_specs <- list()
  design_full <- make_design(meta)
  model_specs[["baseline_full"]] <- list(meta = meta, expr = expr, design = design_full,
    diag = design_diagnostics(design_full, meta, "baseline_full", "baseline", "none", 0L, "full_31v9"))
  for (ref in names(fractions)) {
    pc_obj <- composition_pcs(fractions[[ref]]$fraction, 5L)
    for (k in 2:5) {
      model_id <- paste0(tolower(ref), "_PC1_", k, "_full")
      design <- make_design(meta, pc_obj$scores[, seq_len(k), drop = FALSE])
      model_specs[[model_id]] <- list(meta = meta, expr = expr, design = design,
        diag = design_diagnostics(design, meta, model_id, "composition", ref, k, "full_31v9",
          baseline_vif = NA_real_, original_main = ref == "MuSiC" && k == 5L))
      model_specs[[model_id]]$pc_scores <- pc_obj$scores[, seq_len(k), drop = FALSE]
    }
  }
  baseline_vif_full <- model_specs[["baseline_full"]]$diag$group_VIF[1]
  for (nm in names(model_specs)) if (nm != "baseline_full") {
    model_specs[[nm]]$diag$baseline_group_VIF <- baseline_vif_full
    model_specs[[nm]]$diag$group_VIF_ratio_vs_baseline <-
      model_specs[[nm]]$diag$group_VIF / baseline_vif_full
  }
  tissue_base <- make_design(meta_tissue)
  model_specs[["baseline_tissue"]] <- list(meta = meta_tissue, expr = expr_tissue, design = tissue_base,
    diag = design_diagnostics(tissue_base, meta_tissue, "baseline_tissue", "baseline", "none", 0L, "tissue_31v8",
      baseline_vif = NA_real_))
  tissue_frac <- fractions$MuSiC$fraction[meta_tissue$sample, , drop = FALSE]
  tissue_pc <- composition_pcs(tissue_frac, 5L)
  for (k in 2:5) {
    model_id <- paste0("music_PC1_", k, "_tissue")
    design <- make_design(meta_tissue, tissue_pc$scores[, seq_len(k), drop = FALSE])
    model_specs[[model_id]] <- list(meta = meta_tissue, expr = expr_tissue, design = design,
      diag = design_diagnostics(design, meta_tissue, model_id, "composition", "MuSiC", k, "tissue_31v8",
        baseline_vif = model_specs[["baseline_tissue"]]$diag$group_VIF[1], original_main = k == 5L))
  }
  diag_tab <- do.call(rbind, lapply(model_specs, function(x) x$diag))
  write_csv(diag_tab, file.path(ANALYSIS, "bulk_design_diagnostics.csv"))
  log_note("bulk-model designs assembled: ", nrow(diag_tab), " models; accepted by numeric gate: ", sum(diag_tab$accepted_by_original_numeric_gate))

  # Exact equivalence check between old ~0+group+batch / Tumor-Control and new ~group+batch.
  old_dat <- data.frame(group = factor(as.character(meta$group), levels = c("Control", "Tumor")), batch = meta$batch)
  old_design <- stats::model.matrix(~ 0 + group + batch, data = old_dat)
  old_fit <- limma::lmFit(expr, old_design)
  old_contrast <- limma::makeContrasts(Tumor_minus_Control = groupTumor - groupControl, levels = old_design)
  old_fit <- limma::eBayes(limma::contrasts.fit(old_fit, old_contrast), trend = TRUE, robust = TRUE)
  new_fit <- fit_limma(expr, meta, design_full, model_specs[["baseline_full"]]$diag, "baseline_full")$fit
  old_beta <- old_fit$coefficients[, "Tumor_minus_Control"]
  new_beta <- new_fit$coefficients[, "groupTumor"]
  eq <- data.frame(check = c("coefficient max abs difference", "ordinary residual sigma max abs difference", "raw SE max abs difference", "moderated SE max abs difference", "moderated t max abs difference", "P value max abs difference"),
    value = c(max(abs(old_beta - new_beta)), max(abs(old_fit$sigma - new_fit$sigma)),
      max(abs(old_fit$sigma * old_fit$stdev.unscaled[, 1] - new_fit$sigma * new_fit$stdev.unscaled[, "groupTumor"])),
      max(abs(sqrt(old_fit$s2.post) * old_fit$stdev.unscaled[, 1] - sqrt(new_fit$s2.post) * new_fit$stdev.unscaled[, "groupTumor"])),
      max(abs(old_fit$t[, 1] - new_fit$t[, "groupTumor"])), max(abs(old_fit$p.value[, 1] - new_fit$p.value[, "groupTumor"]))),
    stringsAsFactors = FALSE)
  write_csv(eq, file.path(ANALYSIS, "bulk_naive_design_equivalence.csv"))
  assert_true(max(eq$value) < 1e-10, "old group contrast and new group coefficient are algebraically equivalent")

  log_note("Fitting fresh limma gene tables for all full-rank bulk-model designs")
  fit_tables <- list(); fit_objects <- list()
  for (model_id in names(model_specs)) {
    spec <- model_specs[[model_id]]
    if (!isTRUE(spec$diag$full_rank[1])) {
      log_note("Diagnostic design rank deficient; no fit emitted: ", model_id)
      next
    }
    fitted <- fit_limma(spec$expr, spec$meta, spec$design, spec$diag, model_id)
    if (is.null(fitted)) next
    fit_tables[[model_id]] <- fitted$table; fit_objects[[model_id]] <- fitted$fit
    if (model_id %in% c("baseline_full", "music_PC1_5_full", "baseline_tissue", "music_PC1_5_tissue")) {
      write_csv(fitted$table, file.path(ANALYSIS, paste0("bulk_gene_level_", model_id, ".csv.gz")))
      log_note("Wrote early gene table: ", model_id, " n=", nrow(fitted$table))
    }
  }
  assert_true("baseline_full" %in% names(fit_tables), "fresh full baseline fit completed")
  assert_true("music_PC1_5_full" %in% names(fit_tables), "fresh main MuSiC 5PC fit completed")
  base <- fit_tables$baseline_full
  base_sig <- !is.na(base$adj.P.Val) & base$adj.P.Val < 0.05 & abs(base$logFC) > 1
  summary_rows <- list(); threshold_rows <- list()
  for (model_id in names(fit_tables)) {
    spec <- model_specs[[model_id]]
    if (model_id %in% c("baseline_full", "baseline_tissue")) next
    if (spec$diag$analysis_set[1] == "full_31v9") {
      summary_rows[[model_id]] <- summarize_model(spec$diag, fit_tables[[model_id]], base, base_sig)
      threshold_rows[[model_id]] <- make_threshold_status(fit_tables[[model_id]], base, model_id, "full_31v9")
    }
  }
  summary_tab <- do.call(rbind, summary_rows)
  write_csv(summary_tab, file.path(ANALYSIS, "bulk_model_summary_full.csv"))
  write_csv(do.call(rbind, threshold_rows), file.path(ANALYSIS, "bulk_baseline_threshold_status_by_model.csv"))
  write_csv(do.call(rbind, fit_tables), file.path(ANALYSIS, "bulk_gene_level_all_fits.csv.gz"))
  tissue_base_sig <- !is.na(fit_tables$baseline_tissue$adj.P.Val) & fit_tables$baseline_tissue$adj.P.Val < 0.05 & abs(fit_tables$baseline_tissue$logFC) > 1
  tissue_summaries <- list(); tissue_thresholds <- list()
  for (model_id in names(fit_tables)[grepl("_tissue$", names(fit_tables))]) {
    tissue_summaries[[model_id]] <- summarize_model(model_specs[[model_id]]$diag, fit_tables[[model_id]], fit_tables$baseline_tissue, tissue_base_sig)
    tissue_thresholds[[model_id]] <- make_threshold_status(fit_tables[[model_id]], fit_tables$baseline_tissue, model_id, "tissue_31v8")
  }
  write_csv(do.call(rbind, tissue_summaries), file.path(ANALYSIS, "bulk_tissue_only_model_summary.csv"))
  write_csv(do.call(rbind, tissue_thresholds), file.path(ANALYSIS, "bulk_tissue_only_threshold_status.csv"))
  full_deg_genes <- base$gene[base_sig]
  tissue_map <- match(full_deg_genes, fit_tables$baseline_tissue$gene)
  tissue_compare <- data.frame(comparison = "fresh tissue-only naive vs fresh full baseline; fixed fresh full baseline DEG universe",
    full_baseline_DEGs = length(full_deg_genes),
    retained_tissue_thresholds = sum(fit_tables$baseline_tissue$adj.P.Val[tissue_map] < 0.05 & abs(fit_tables$baseline_tissue$logFC[tissue_map]) > 1, na.rm = TRUE),
    stringsAsFactors = FALSE)
  write_csv(tissue_compare, file.path(ANALYSIS, "bulk_tissue_only_vs_full_baseline.csv"))

  # Group composition overlap and sample-level plotting/inspection source.
  batch_tab <- as.data.frame(table(batch = meta$batch, group = meta$group), stringsAsFactors = FALSE)
  write_csv(batch_tab, file.path(ANALYSIS, "bulk_batch_by_group.csv"))
  overlap_rows <- list(); distribution_rows_all <- list(); score_rows <- list()
  for (ref in names(fractions)) {
    f <- fractions[[ref]]$fraction
    pc <- composition_pcs(f, 5L)$scores
    sc <- data.frame(sample = meta$sample, group = meta$group, batch = meta$batch, reference = ref,
      as.data.frame(f, check.names = FALSE), as.data.frame(pc, check.names = FALSE), check.names = FALSE)
    score_rows[[ref]] <- sc
    for (feature in c(colnames(f), colnames(pc))) {
      val <- if (feature %in% colnames(f)) f[, feature] else pc[, feature]
      names(val) <- rownames(f)
      distribution_rows_all[[paste(ref, feature)]] <- distribution_rows(val[meta$sample], meta, feature, ref)
      overlap_rows[[paste(ref, feature)]] <- range_overlap(val[meta$group == "Control"], val[meta$group == "Tumor"], feature, "group", ref)
    }
  }
  write_csv(do.call(rbind, score_rows), file.path(ANALYSIS, "bulk_sample_composition_scores.csv"))
  write_csv(do.call(rbind, distribution_rows_all), file.path(ANALYSIS, "bulk_composition_distribution.csv"))
  write_csv(do.call(rbind, overlap_rows), file.path(ANALYSIS, "bulk_composition_range_overlap.csv"))
  log_note("bulk-model complete: fresh baseline DEG n=", sum(base_sig), "; fresh MuSiC PC1-5 DEG n=", sum(fit_tables$music_PC1_5_full$adj.P.Val < 0.05 & abs(fit_tables$music_PC1_5_full$logFC) > 1, na.rm = TRUE))

  # Historical tables are admitted only after all fresh fits have completed.
  historical_path <- file.path(ROOT, "results/deconv/deg_composition_adjusted.csv")
  historical <- read_csv(historical_path)
  hist_idx <- match(rownames(expr), historical$gene)
  assert_true(!anyNA(hist_idx) && nrow(historical) == nrow(expr), "historical comparison covers full discovery gene universe")
  fresh_base <- fit_tables$baseline_full
  fresh_music <- fit_tables$music_PC1_5_full
  hist_gene <- data.frame(gene = rownames(expr),
    historical_original_logFC = historical$original_logFC[hist_idx],
    historical_original_adjP = historical$original_adj.P.Val[hist_idx],
    historical_original_sig = as.logical(historical$original_sig[hist_idx]),
    fresh_group_batch_logFC = fresh_base$logFC,
    fresh_group_batch_adjP = fresh_base$adj.P.Val,
    fresh_group_batch_sig = !is.na(fresh_base$adj.P.Val) & fresh_base$adj.P.Val < .05 & abs(fresh_base$logFC) > 1,
    historical_adjusted_logFC = historical$adjusted_logFC[hist_idx],
    historical_adjusted_adjP = historical$adjusted_adj.P.Val[hist_idx],
    historical_adjusted_sig = as.logical(historical$adjusted_sig[hist_idx]),
    fresh_MuSiC_PC1_5_logFC = fresh_music$logFC,
    fresh_MuSiC_PC1_5_adjP = fresh_music$adj.P.Val,
    fresh_MuSiC_PC1_5_sig = !is.na(fresh_music$adj.P.Val) & fresh_music$adj.P.Val < .05 & abs(fresh_music$logFC) > 1,
    stringsAsFactors = FALSE)
  write_csv(hist_gene, file.path(ANALYSIS, "bulk_postfit_historical_gene_comparison.csv.gz"))
  historical_summary <- do.call(rbind, lapply(list(
    historical_original_vs_fresh_group_batch = c("historical_original_logFC", "fresh_group_batch_logFC",
      "historical_original_sig", "fresh_group_batch_sig"),
    historical_adjusted_vs_fresh_MuSiC_PC1_5 = c("historical_adjusted_logFC", "fresh_MuSiC_PC1_5_logFC",
      "historical_adjusted_sig", "fresh_MuSiC_PC1_5_sig")
  ), function(columns) {
    nm <- if (columns[1] == "historical_original_logFC") "historical_original_vs_fresh_group_batch" else "historical_adjusted_vs_fresh_MuSiC_PC1_5"
    x <- hist_gene[[columns[1]]]; y <- hist_gene[[columns[2]]]
    old_flag <- hist_gene[[columns[3]]]; new_flag <- hist_gene[[columns[4]]]
    met <- effect_metrics(x, y)
    data.frame(comparison = nm, genes_compared = sum(is.finite(x) & is.finite(y)),
      max_abs_effect_difference = max(abs(x - y), na.rm = TRUE),
      historical_DEG_count = sum(old_flag, na.rm = TRUE), fresh_DEG_count = sum(new_flag, na.rm = TRUE),
      membership_discordant = sum(old_flag != new_flag, na.rm = TRUE),
      historical_only = sum(old_flag & !new_flag, na.rm = TRUE), fresh_only = sum(!old_flag & new_flag, na.rm = TRUE),
      Pearson = met[["Pearson"]], Spearman = met[["Spearman"]], CCC = met[["CCC"]],
      stringsAsFactors = FALSE)
  }))
  write_csv(historical_summary, file.path(ANALYSIS, "bulk_postfit_historical_comparison_summary.csv"))
  write_csv(data.frame(path = historical_path, sha256 = sha256_file(historical_path),
    role = "read only after independent fresh fits; postfit comparison only", stringsAsFactors = FALSE),
    file.path(ANALYSIS, "bulk_postfit_historical_input_hash.csv"))
  log_note("Postfit historical comparison completed only after fresh bulk-model tables; historical source remains read-only")

  list(paths = paths, expr = expr, meta = meta, ext_expr = ext_expr, ext_meta = ext_meta, hubs = hubs,
    hubtab = hubtab, fit_tables = fit_tables, model_specs = model_specs, diag = diag_tab,
    fractions = fractions, input_hash = input_hash, new_fit = new_fit, old_fit = old_fit)
}

run_agreement <- function(state) {
  log_note("effect-agreement: validating subtype evidence and fitting independent external group-only comparisons")
  ext_meta <- state$ext_meta
  subtype_evidence <- read_csv(state$paths$subtype_evidence)
  ev_idx <- match(ext_meta$sample, subtype_evidence$sample)
  assert_true(!anyNA(ev_idx), "subtype evidence covers every GSE108524 sample")
  ext_meta$assigned_subtype <- as.character(subtype_evidence$assigned_subtype[ev_idx])
  ext_meta$source_disease_state <- as.character(ext_meta$disease_state)
  ext_meta$evidence_title <- as.character(subtype_evidence$title[ev_idx])
  ext_meta$classification_basis <- as.character(subtype_evidence$classification_basis[ev_idx])
  subtype_counts <- table(ext_meta$assigned_subtype)
  assert_true(identical(as.integer(subtype_counts[c("control", "sporadic", "NF2")]), c(4L, 10L, 17L)),
    paste("expected evidence subtype counts control/sporadic/NF2=4/10/17; got",
      paste(names(subtype_counts), as.integer(subtype_counts), collapse = "; ")))
  assert_true(all(ext_meta$assigned_subtype[ext_meta$group == "Control"] == "control") &&
    all(ext_meta$assigned_subtype[ext_meta$group == "Tumor"] %in% c("sporadic", "NF2")),
    "subtype evidence and processed group labels align")
  write_csv(ext_meta[, c("sample", "group", "disease_state", "source_disease_state", "assigned_subtype",
    "evidence_title", "classification_basis")], file.path(ANALYSIS, "agreement_subtype_evidence_alignment.csv"))
  log_note("Subtype evidence confirms 4 controls, 10 sporadic, 17 NF2; sporadic raw disease_state label retained verbatim")

  genes_common <- rownames(state$expr)[rownames(state$expr) %in% rownames(state$ext_expr)]
  assert_true(length(genes_common) >= length(state$hubs) && all(state$hubs %in% genes_common), "common measurable gene universe covers fixed 346")
  discovery_models <- list(
    fresh_group_batch = state$fit_tables$baseline_full,
    fresh_MuSiC_PC1_5 = state$fit_tables$music_PC1_5_full
  )
  external_defs <- list(
    all_VS = ext_meta$sample[ext_meta$group == "Tumor"],
    sporadic = ext_meta$sample[ext_meta$assigned_subtype == "sporadic"],
    NF2 = ext_meta$sample[ext_meta$assigned_subtype == "NF2"]
  )
  control_ids <- ext_meta$sample[ext_meta$group == "Control"]
  external_fits <- list()
  for (contrast_id in names(external_defs)) {
    external_fits[[contrast_id]] <- external_group_fit(state$ext_expr[genes_common, , drop = FALSE], ext_meta,
      contrast_id, external_defs[[contrast_id]], control_ids)
  }
  write_csv(do.call(rbind, lapply(external_fits, function(x) x$table)),
    file.path(ANALYSIS, "agreement_external_gene_level_group_only.csv.gz"))
  common_tab <- list(); agreement_rows <- list(); fdr_rows <- list()
  for (contrast_id in names(external_fits)) {
    ext_tab <- external_fits[[contrast_id]]$table
    ext_vec <- stats::setNames(ext_tab$logFC, ext_tab$gene)
    ext_fdr <- stats::setNames(ext_tab$adj.P.Val, ext_tab$gene)
    for (model_name in names(discovery_models)) {
      disc_tab <- discovery_models[[model_name]]
      disc_vec <- stats::setNames(disc_tab$logFC, disc_tab$gene)
      disc_fdr <- stats::setNames(disc_tab$adj.P.Val, disc_tab$gene)
      for (scope in c("fixed_submitted_346", "all_common_measurable")) {
        use_genes <- if (scope == "fixed_submitted_346") state$hubs else genes_common
        m <- effect_metrics(disc_vec[use_genes], ext_vec[use_genes])
        row <- data.frame(comparison_id = paste0(contrast_id, "_vs_", model_name),
          external_contrast = contrast_id, discovery_effect_model = model_name, scope = scope,
          n_universe = length(use_genes), as.list(m), stringsAsFactors = FALSE, check.names = FALSE)
        agreement_rows[[paste(contrast_id, model_name, scope)]] <- row
        fdr_rows[[paste(contrast_id, model_name, scope)]] <- fdr_rank_stat(use_genes, disc_fdr, ext_fdr,
          paste0(contrast_id, "_vs_", model_name), scope)
      }
    }
  }
  for (model_name in names(discovery_models)) {
    a <- discovery_models$fresh_group_batch
    b <- discovery_models$fresh_MuSiC_PC1_5
    va <- stats::setNames(a$logFC, a$gene)
    vb <- stats::setNames(b$logFC, b$gene)
    fa <- stats::setNames(a$adj.P.Val, a$gene)
    fb <- stats::setNames(b$adj.P.Val, b$gene)
    for (scope in c("fixed_submitted_346", "all_common_measurable")) {
      use_genes <- if (scope == "fixed_submitted_346") state$hubs else genes_common
      metrics <- effect_metrics(va[use_genes], vb[use_genes])
      agreement_rows[[paste("discovery", scope)]] <- data.frame(comparison_id = "fresh_group_batch_vs_fresh_MuSiC_PC1_5",
        external_contrast = "within_discovery", discovery_effect_model = "fresh_naive_vs_fresh_adjusted",
        scope = scope, n_universe = length(use_genes), as.list(metrics), stringsAsFactors = FALSE, check.names = FALSE)
      fdr_rows[[paste("discovery", scope)]] <- fdr_rank_stat(use_genes, fa, fb,
        "fresh_group_batch_vs_fresh_MuSiC_PC1_5", scope)
    }
    break
  }
  point_tab <- do.call(rbind, agreement_rows)
  write_csv(point_tab, file.path(ANALYSIS, "agreement_effect_agreement_point_estimates.csv"))
  write_csv(do.call(rbind, fdr_rows), file.path(ANALYSIS, "agreement_FDR_rank_agreement_separate.csv"))
  log_note("effect-agreement point estimates complete: matched all-common gene universe n=", length(genes_common),
    "; fixed primary denominator n=", length(state$hubs))

  log_note("Starting fixed 500-draw stratified sample bootstrap; seed=", SEED,
    "; external control draw shared among all three comparisons within each replicate")
  set.seed(SEED)
  discovery_strata_raw <- split(seq_len(nrow(state$meta)),
    interaction(state$meta$group, state$meta$batch, drop = TRUE, sep = "__"))
  discovery_strata <- lapply(discovery_strata_raw, function(ix) { z <- rep(FALSE, nrow(state$meta)); z[ix] <- TRUE; z })
  external_strata <- list(
    Control = ext_meta$assigned_subtype == "control",
    sporadic = ext_meta$assigned_subtype == "sporadic",
    NF2 = ext_meta$assigned_subtype == "NF2"
  )
  bootstrap_rows <- vector("list", B * 7L)
  draw_rows <- vector("list", B * 3L)
  record_pos <- 0L; draw_pos <- 0L
  bootstrap_diag_rows <- vector("list", B)
  for (replicate in seq_len(B)) {
    dd <- sample_stratified_draw(state$meta, discovery_strata)
    ed <- sample_stratified_draw(ext_meta, external_strata)
    d_ids_orig <- state$meta$sample[dd$index]
    d_boot_ids <- sprintf("D%03d_%03d", replicate, seq_along(dd$index))
    dmeta <- state$meta[dd$index, , drop = FALSE]; dmeta$sample <- d_boot_ids
    dexpr <- state$expr[, d_ids_orig, drop = FALSE]; colnames(dexpr) <- d_boot_ids
    dfrac <- state$fractions$MuSiC$fraction[d_ids_orig, , drop = FALSE]; rownames(dfrac) <- d_boot_ids
    ddesign_naive <- make_design(dmeta)
    ddiag_naive <- design_diagnostics(ddesign_naive, dmeta, paste0("boot", replicate, "_naive"),
      "bootstrap", "none", 0L, "bootstrap_discovery")
    dcoef_naive <- matrix_ols_group_coef(dexpr, dmeta$sample, ddesign_naive)
    ddesign_adjusted <- NULL; ddiag_adjusted <- NULL; dcoef_adjusted <- NULL
    tryCatch({
      dpc <- composition_pcs(dfrac, 5L)$scores
      ddesign_adjusted <- make_design(dmeta, dpc)
      ddiag_adjusted <- design_diagnostics(ddesign_adjusted, dmeta, paste0("boot", replicate, "_MuSiC_PC1_5"),
        "bootstrap", "MuSiC", 5L, "bootstrap_discovery", original_main = TRUE)
      dcoef_adjusted <- matrix_ols_group_coef(dexpr, dmeta$sample, ddesign_adjusted)
    }, error = function(e) {
      ddiag_adjusted <<- data.frame(full_rank = FALSE, accepted_by_original_numeric_gate = FALSE,
        gate_failure_reason = paste0("PCA/design error: ", conditionMessage(e)), rank = 0L,
        n_columns = 8L, residual_df = NA_integer_, condition_number = Inf,
        original_approx_condition_number = Inf, exact_2norm_condition_number = Inf,
        group_VIF = Inf, group_R2_from_nuisance = NA_real_)
    })
    control_idx <- ed$by_stratum$Control
    sporadic_idx <- ed$by_stratum$sporadic
    nf2_idx <- ed$by_stratum$NF2
    ext_sample_effect <- function(tumor_idx) {
      tumor_ids <- ext_meta$sample[tumor_idx]
      control_ids_draw <- ext_meta$sample[control_idx]
      as.numeric(rowMeans(state$ext_expr[state$hubs, tumor_ids, drop = FALSE]) -
        rowMeans(state$ext_expr[state$hubs, control_ids_draw, drop = FALSE]))
    }
    names_all <- state$ext_meta$sample[ed$by_stratum$sporadic]
    names_nf2 <- state$ext_meta$sample[ed$by_stratum$NF2]
    all_tumor_idx <- c(sporadic_idx, nf2_idx)
    ext_effects <- list(all_VS = ext_sample_effect(all_tumor_idx),
      sporadic = ext_sample_effect(sporadic_idx), NF2 = ext_sample_effect(nf2_idx))
    names(ext_effects$all_VS) <- state$hubs; names(ext_effects$sporadic) <- state$hubs; names(ext_effects$NF2) <- state$hubs
    discovery_coefficients <- list(fresh_group_batch = dcoef_naive, fresh_MuSiC_PC1_5 = dcoef_adjusted)
    for (contrast_id in names(ext_effects)) {
      for (model_name in names(discovery_coefficients)) {
        x <- discovery_coefficients[[model_name]]
        y <- ext_effects[[contrast_id]]
        x_hub <- if (!is.null(x) && !is.null(names(x)) && all(state$hubs %in% names(x))) x[state$hubs] else numeric()
        estimable <- length(x_hub) == length(state$hubs) && all(is.finite(x_hub)) && all(is.finite(y))
        comparison_id <- paste0(contrast_id, "_vs_", model_name)
        gate_failure <- if (model_name == "fresh_MuSiC_PC1_5") {
          is.null(ddiag_adjusted) || !isTRUE(ddiag_adjusted$accepted_by_original_numeric_gate[1])
        } else !isTRUE(ddiag_naive$accepted_by_original_numeric_gate[1])
        met <- if (estimable) effect_metrics(x_hub, y[state$hubs]) else numeric()
        reason <- if (estimable) NULL else if (model_name == "fresh_MuSiC_PC1_5" && (is.null(ddiag_adjusted) || !isTRUE(ddiag_adjusted$full_rank[1]))) "composition fit unestimable/rank deficient" else "OLS coefficient vector unestimable"
        record_pos <- record_pos + 1L
        bootstrap_rows[[record_pos]] <- emit_boot_metric_rows(replicate, comparison_id, met,
          gate_failure, estimable, reason)
      }
    }
    naive_hub <- if (!is.null(dcoef_naive) && !is.null(names(dcoef_naive)) && all(state$hubs %in% names(dcoef_naive))) dcoef_naive[state$hubs] else numeric()
    adjusted_hub <- if (!is.null(dcoef_adjusted) && !is.null(names(dcoef_adjusted)) && all(state$hubs %in% names(dcoef_adjusted))) dcoef_adjusted[state$hubs] else numeric()
    disc_estimable <- length(naive_hub) == length(state$hubs) && length(adjusted_hub) == length(state$hubs) &&
      all(is.finite(naive_hub)) && all(is.finite(adjusted_hub))
    disc_met <- if (disc_estimable) effect_metrics(naive_hub, adjusted_hub) else numeric()
    record_pos <- record_pos + 1L
    bootstrap_rows[[record_pos]] <- emit_boot_metric_rows(replicate, "fresh_group_batch_vs_fresh_MuSiC_PC1_5",
      disc_met, is.null(ddiag_adjusted) || !isTRUE(ddiag_adjusted$accepted_by_original_numeric_gate[1]),
      disc_estimable, if (disc_estimable) NULL else "one or both discovery OLS coefficient vectors unestimable")

    # Auditable draw-index table: duplicates are retained as distinct draw positions.
    for (nm in names(dd$by_stratum)) {
      ids <- dd$by_stratum[[nm]]
      draw_pos <- draw_pos + 1L
      draw_rows[[draw_pos]] <- data.frame(replicate = replicate, cohort = "GSE39645", stratum = nm,
        draw_position = seq_along(ids), source_index = ids, sample = state$meta$sample[ids], stringsAsFactors = FALSE)
    }
    for (nm in names(ed$by_stratum)) {
      ids <- ed$by_stratum[[nm]]
      draw_pos <- draw_pos + 1L
      draw_rows[[draw_pos]] <- data.frame(replicate = replicate, cohort = "GSE108524", stratum = nm,
        draw_position = seq_along(ids), source_index = ids, sample = ext_meta$sample[ids], stringsAsFactors = FALSE)
    }
    naive_vif_fail <- !is.finite(ddiag_naive$group_VIF[1]) || ddiag_naive$group_VIF[1] > VIF_LIMIT
    naive_approx_fail <- !is.finite(ddiag_naive$original_approx_condition_number[1]) ||
      ddiag_naive$original_approx_condition_number[1] > CONDITION_LIMIT
    naive_rank_fail <- !isTRUE(ddiag_naive$full_rank[1])
    adjusted_rank_fail <- is.null(ddiag_adjusted) || !isTRUE(ddiag_adjusted$full_rank[1])
    adjusted_vif <- if (is.null(ddiag_adjusted)) NA_real_ else ddiag_adjusted$group_VIF[1]
    adjusted_approx <- if (is.null(ddiag_adjusted)) NA_real_ else ddiag_adjusted$original_approx_condition_number[1]
    adjusted_exact <- if (is.null(ddiag_adjusted)) NA_real_ else ddiag_adjusted$exact_2norm_condition_number[1]
    adjusted_vif_fail <- !is.finite(adjusted_vif) || adjusted_vif > VIF_LIMIT
    adjusted_approx_fail <- !is.finite(adjusted_approx) || adjusted_approx > CONDITION_LIMIT
    bootstrap_diag_rows[[replicate]] <- data.frame(replicate = replicate,
      discovery_naive_rank = ddiag_naive$rank[1], discovery_naive_columns = ddiag_naive$n_columns[1],
      discovery_naive_residual_df = ddiag_naive$residual_df[1],
      discovery_naive_group_VIF = ddiag_naive$group_VIF[1],
      discovery_naive_original_approx_condition = ddiag_naive$original_approx_condition_number[1],
      discovery_naive_exact_2norm_condition = ddiag_naive$exact_2norm_condition_number[1],
      discovery_naive_rank_failure = naive_rank_fail, discovery_naive_approx_condition_failure = naive_approx_fail,
      discovery_naive_VIF_failure = naive_vif_fail,
      discovery_naive_gate_failure = !isTRUE(ddiag_naive$accepted_by_original_numeric_gate[1]),
      discovery_adjusted_rank = if (is.null(ddiag_adjusted)) NA_integer_ else ddiag_adjusted$rank[1],
      discovery_adjusted_columns = if (is.null(ddiag_adjusted)) NA_integer_ else ddiag_adjusted$n_columns[1],
      discovery_adjusted_residual_df = if (is.null(ddiag_adjusted)) NA_integer_ else ddiag_adjusted$residual_df[1],
      discovery_adjusted_group_VIF = adjusted_vif,
      discovery_adjusted_original_approx_condition = adjusted_approx,
      discovery_adjusted_exact_2norm_condition = adjusted_exact,
      discovery_adjusted_rank_failure = adjusted_rank_fail,
      discovery_adjusted_approx_condition_failure = adjusted_approx_fail,
      discovery_adjusted_VIF_failure = adjusted_vif_fail,
      discovery_adjusted_gate_failure = is.null(ddiag_adjusted) || !isTRUE(ddiag_adjusted$accepted_by_original_numeric_gate[1]),
      discovery_adjusted_gate_reason = if (is.null(ddiag_adjusted)) "design unavailable" else ddiag_adjusted$gate_failure_reason[1],
      discovery_adjusted_estimable = !is.null(dcoef_adjusted), shared_control_n = length(control_idx),
      sporadic_n = length(sporadic_idx), NF2_n = length(nf2_idx), all_VS_tumor_n = length(all_tumor_idx),
      stringsAsFactors = FALSE)
    if (replicate %% 50L == 0L) log_note("Bootstrap completed ", replicate, "/", B,
      "; metric rows recorded=", record_pos, "; draw-index batches=", draw_pos)
  }
  bootstrap_table <- do.call(rbind, bootstrap_rows[seq_len(record_pos)])
  bootstrap_draw_table <- do.call(rbind, draw_rows[seq_len(draw_pos)])
  bootstrap_diag_table <- do.call(rbind, bootstrap_diag_rows)
  write_csv(bootstrap_table, file.path(ANALYSIS, "agreement_bootstrap_metric_draws.csv.gz"))
  write_csv(bootstrap_draw_table, file.path(ANALYSIS, "agreement_bootstrap_resample_indices.csv.gz"))
  write_csv(bootstrap_diag_table, file.path(ANALYSIS, "agreement_bootstrap_design_diagnostics.csv"))
  gate_counts <- data.frame(
    failure_type = c("rank deficient", "original approximate kappa > 1000/nonfinite", "group VIF > 10/nonfinite", "any original numeric gate failure"),
    discovery_draws = c(sum(bootstrap_diag_table$discovery_adjusted_rank_failure),
      sum(bootstrap_diag_table$discovery_adjusted_approx_condition_failure),
      sum(bootstrap_diag_table$discovery_adjusted_VIF_failure),
      sum(bootstrap_diag_table$discovery_adjusted_gate_failure)),
    requested_B = B, stringsAsFactors = FALSE)
  write_csv(gate_counts, file.path(ANALYSIS, "agreement_bootstrap_gate_failure_counts.csv"))
  gate_patterns <- as.data.frame(table(
    rank_failure = bootstrap_diag_table$discovery_adjusted_rank_failure,
    approximate_condition_failure = bootstrap_diag_table$discovery_adjusted_approx_condition_failure,
    VIF_failure = bootstrap_diag_table$discovery_adjusted_VIF_failure), stringsAsFactors = FALSE)
  write_csv(gate_patterns, file.path(ANALYSIS, "agreement_bootstrap_gate_failure_patterns.csv"))
  boot_summary <- summarize_bootstrap(bootstrap_table, B, SEED)
  point_primary <- point_tab[point_tab$scope == "fixed_submitted_346", , drop = FALSE]
  point_metric_names <- unique(bootstrap_table$metric)
  point_long <- do.call(rbind, lapply(point_metric_names, function(metric_name) {
    if (!metric_name %in% names(point_primary)) return(NULL)
    data.frame(comparison_id = point_primary$comparison_id, metric = metric_name,
      point_estimate = as.numeric(point_primary[[metric_name]]), stringsAsFactors = FALSE)
  }))
  boot_summary$point_estimate <- point_long$point_estimate[
    match(paste(boot_summary$comparison_id, boot_summary$metric),
      paste(point_long$comparison_id, point_long$metric))]
  write_csv(boot_summary, file.path(ANALYSIS, "agreement_bootstrap_percentile_intervals.csv"))
  log_note("effect-agreement bootstrap complete: finite rows=", sum(bootstrap_table$valid), "/", nrow(bootstrap_table),
    "; adjusted discovery gate failures=", sum(bootstrap_diag_table$discovery_adjusted_gate_failure))
  list(external_fits = external_fits, point = point_tab, fdr = do.call(rbind, fdr_rows),
    bootstrap = bootstrap_table, bootstrap_summary = boot_summary, bootstrap_draws = bootstrap_draw_table,
    bootstrap_diagnostics = bootstrap_diag_table, subtype_meta = ext_meta, common_genes = genes_common)
}

make_figures <- function(bulk_state, agreement_state) {
  log_note("Creating S15/S16 from fresh CSVs; Arial and approved pastel palette")
  palette_ref <- c(MuSiC = "#e7c5c4", NNLS = "#85bcbc", `GSE216783 NNLS` = "#9fb3d8")
  palette_model <- c(fresh_group_batch = "#e7c5c4", fresh_MuSiC_PC1_5 = "#85bcbc")
  s15_data <- read_csv(file.path(ANALYSIS, "bulk_model_summary_full.csv"))
  s15_diag <- read_csv(file.path(ANALYSIS, "bulk_design_diagnostics.csv"))
  s15_data$retained_percent <- 100 * s15_data$both_thresholds_retained / s15_data$baseline_fresh_DEG_n
  s15_data$reference_label <- ifelse(s15_data$reference == "GSE216783_reference", "GSE216783 NNLS", s15_data$reference)
  s15_data$label <- paste0(s15_data$reference_label, " PC1-", s15_data$pc_n)
  s15_data$reference_label <- factor(s15_data$reference_label, levels = names(palette_ref))
  s15_data$label <- factor(s15_data$label, levels = rev(unique(s15_data$label)))
  write_csv(s15_data, file.path(ANALYSIS, "Figure_S15_effect_threshold_source.csv"))
  p15a <- ggplot2::ggplot(s15_data, ggplot2::aes(x = retained_percent, y = label, color = reference_label)) +
    ggplot2::geom_vline(xintercept = 100, color = "#555555", linetype = "dashed", linewidth = .35) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(values = palette_ref, drop = FALSE) +
    ggplot2::labs(x = "Baseline DEGs retained at both thresholds (%)", y = NULL, color = "Reference") +
    theme_publication() + ggplot2::theme(legend.position = "bottom", axis.text.y = ggplot2::element_text(size = 5.8))

  p15b <- ggplot2::ggplot(s15_data, ggplot2::aes(x = median_absolute_logFC_change_baseline_DEGs, y = label, color = reference_label)) +
    ggplot2::geom_vline(xintercept = 0, color = "#555555", linewidth = .3) +
    ggplot2::geom_point(size = 2) + ggplot2::scale_color_manual(values = palette_ref, drop = FALSE) +
    ggplot2::labs(x = "Median |Δ log2FC| among baseline DEGs", y = NULL, color = "Reference") +
    theme_publication() + ggplot2::theme(legend.position = "none", axis.text.y = ggplot2::element_text(size = 5.8))

  se_long <- rbind(
    data.frame(model_id = s15_data$model_id, label = s15_data$label, reference_label = s15_data$reference_label,
      estimate = s15_data$median_raw_SE_ratio_baseline_DEGs, precision = "Raw SE ratio", stringsAsFactors = FALSE),
    data.frame(model_id = s15_data$model_id, label = s15_data$label, reference_label = s15_data$reference_label,
      estimate = s15_data$median_moderated_SE_ratio_baseline_DEGs, precision = "Moderated SE ratio", stringsAsFactors = FALSE)
  )
  se_long$precision <- factor(se_long$precision, levels = c("Raw SE ratio", "Moderated SE ratio"))
  write_csv(se_long, file.path(ANALYSIS, "Figure_S15_SE_ratio_source.csv"))
  p15c <- ggplot2::ggplot(se_long, ggplot2::aes(x = estimate, y = label, color = reference_label, shape = precision)) +
    ggplot2::geom_vline(xintercept = 1, color = "#555555", linetype = "dashed", linewidth = .3) +
    ggplot2::geom_point(position = ggplot2::position_dodge(width = .45), size = 1.7) +
    ggplot2::scale_color_manual(values = palette_ref, drop = FALSE) +
    ggplot2::labs(x = "Median gene-wise SE ratio vs baseline", y = NULL, color = "Reference", shape = "SE") +
    theme_publication() + ggplot2::theme(legend.position = "bottom", axis.text.y = ggplot2::element_text(size = 5.8)) +
    ggplot2::guides(color = "none")

  gate_data <- s15_diag[(s15_diag$analysis_set == "full_31v9" & s15_diag$model_class == "composition") |
    s15_diag$model_id == "music_PC1_5_tissue", , drop = FALSE]
  gate_data$reference_label <- ifelse(gate_data$reference == "GSE216783_reference", "GSE216783 NNLS", gate_data$reference)
  gate_data$label <- ifelse(gate_data$analysis_set == "tissue_31v8",
    "Tissue-only MuSiC PC1-5", paste0(gate_data$reference_label, " PC1-", gate_data$pc_n))
  gate_data$reference_label <- factor(gate_data$reference_label, levels = names(palette_ref))
  gate_data$label <- factor(gate_data$label, levels = c(levels(s15_data$label), "Tissue-only MuSiC PC1-5"))
  write_csv(gate_data, file.path(ANALYSIS, "Figure_S15_gate_source.csv"))
  p15d <- ggplot2::ggplot(gate_data, ggplot2::aes(x = group_VIF, y = label, color = reference_label, shape = accepted_by_original_numeric_gate)) +
    ggplot2::geom_vline(xintercept = VIF_LIMIT, color = "#555555", linetype = "dashed", linewidth = .3) +
    ggplot2::geom_point(size = 2) + ggplot2::scale_color_manual(values = palette_ref, drop = FALSE) +
    ggplot2::scale_shape_manual(values = c("FALSE" = 4, "TRUE" = 16), labels = c("FALSE" = "Fails gate", "TRUE" = "Passes gate")) +
    ggplot2::labs(x = "Group VIF (gate ≤ 10)", y = NULL, color = "Reference", shape = "Original gate") +
    theme_publication() + ggplot2::theme(legend.position = "bottom", axis.text.y = ggplot2::element_text(size = 5.8)) +
    ggplot2::guides(color = "none")
  p15 <- (p15a | p15b) / (p15c | p15d) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(tag_levels = "a", title = "Composition-adjusted discovery sensitivity",
      subtitle = "a–c: 12 full-cohort models (31 VS/9 controls). d: these plus separately labeled tissue-only MuSiC PC5 (31 VS/8 controls; VIF 10.353, fails gate). Gate uses approximate kappa and VIF.",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(family = "Arial", face = "bold", size = 10),
        plot.subtitle = ggplot2::element_text(family = "Arial", size = 7),
        plot.tag = ggplot2::element_text(family = "Arial", face = "bold", size = 10))) &
    ggplot2::theme(legend.position = "bottom")
  write_plot_formats(p15, file.path(FIGURES, "Figure_S15_model_sensitivity"), 11, 8.1)

  naive_tab <- bulk_state$fit_tables$baseline_full
  hub_naive <- naive_tab[match(bulk_state$hubs, naive_tab$gene), c("gene", "logFC"), drop = FALSE]
  names(hub_naive)[2] <- "discovery_logFC"
  ext_tab <- read_csv(file.path(ANALYSIS, "agreement_external_gene_level_group_only.csv.gz"))
  ext_hub <- ext_tab[ext_tab$gene %in% bulk_state$hubs, , drop = FALSE]
  scatter <- merge(ext_hub, hub_naive, by = "gene", all.x = TRUE, sort = FALSE)
  scatter$contrast_label <- factor(scatter$contrast_id, levels = c("all_VS", "sporadic", "NF2"),
    labels = c("All VS (27)", "Sporadic (10)", "NF2-associated (17)"))
  point_stats <- read_csv(file.path(ANALYSIS, "agreement_effect_agreement_point_estimates.csv"))
  point_stats <- point_stats[point_stats$discovery_effect_model == "fresh_group_batch" &
    point_stats$scope == "fixed_submitted_346" & point_stats$external_contrast != "within_discovery", , drop = FALSE]
  point_stats$contrast_label <- factor(point_stats$external_contrast, levels = c("all_VS", "sporadic", "NF2"),
    labels = c("All VS (27)", "Sporadic (10)", "NF2-associated (17)"))
  point_stats$annotation <- sprintf("r = %.2f; direction = %.1f%%", point_stats$Pearson, 100 * point_stats$direction_rate)
  write_csv(scatter, file.path(ANALYSIS, "Figure_S16_scatter_source.csv"))
  scatter_lim <- range(c(scatter$discovery_logFC, scatter$logFC), finite = TRUE)
  p16a <- ggplot2::ggplot(scatter, ggplot2::aes(x = discovery_logFC, y = logFC, color = contrast_label)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, color = "#444444", linetype = "dashed", linewidth = .35) +
    ggplot2::geom_point(alpha = .52, size = .75, stroke = 0) +
    ggplot2::geom_text(data = point_stats, ggplot2::aes(x = -Inf, y = Inf, label = annotation),
      inherit.aes = FALSE, hjust = -.05, vjust = 1.2, size = 2.3, family = "Arial", color = "#333333", show.legend = FALSE) +
    ggplot2::facet_wrap(~contrast_label, nrow = 1) +
    ggplot2::scale_color_manual(values = c("All VS (27)" = "#e7c5c4", "Sporadic (10)" = "#85bcbc", "NF2-associated (17)" = "#cab8e4")) +
    ggplot2::coord_equal(xlim = scatter_lim, ylim = scatter_lim, expand = TRUE) +
    ggplot2::labs(x = "Discovery group+batch OLS log2 effect", y = "External group-only OLS log2 difference") +
    theme_publication() + ggplot2::theme(legend.position = "none", strip.text = ggplot2::element_text(family = "Arial", face = "bold", size = 7))

  ci <- read_csv(file.path(ANALYSIS, "agreement_bootstrap_percentile_intervals.csv"))
  ci_ids <- c("all_VS_vs_fresh_group_batch", "all_VS_vs_fresh_MuSiC_PC1_5",
    "sporadic_vs_fresh_group_batch", "sporadic_vs_fresh_MuSiC_PC1_5",
    "NF2_vs_fresh_group_batch", "NF2_vs_fresh_MuSiC_PC1_5")
  ci <- ci[ci$metric %in% c("direction_rate", "Pearson", "Spearman", "CCC") & ci$comparison_id %in% ci_ids, , drop = FALSE]
  ci$comparison_label <- sub("_vs_fresh_(group_batch|MuSiC_PC1_5)$", "", ci$comparison_id)
  ci$comparison_label <- factor(ci$comparison_label, levels = c("all_VS", "sporadic", "NF2"),
    labels = c("All VS", "Sporadic", "NF2"))
  ci$discovery_model <- ifelse(grepl("MuSiC_PC1_5", ci$comparison_id), "fresh_MuSiC_PC1_5", "fresh_group_batch")
  ci$metric <- factor(ci$metric, levels = c("direction_rate", "Pearson", "Spearman", "CCC"),
    labels = c("Direction match", "Pearson r", "Spearman ρ", "Lin CCC"))
  write_csv(ci, file.path(ANALYSIS, "Figure_S16_interval_source.csv"))
  p16b <- ggplot2::ggplot(ci, ggplot2::aes(x = comparison_label, y = point_estimate, ymin = percentile_2.5,
    ymax = percentile_97.5, color = discovery_model, group = discovery_model)) +
    ggplot2::geom_hline(yintercept = 0, color = "#666666", linetype = "dashed", linewidth = .25) +
    ggplot2::geom_pointrange(position = ggplot2::position_dodge(width = .5), linewidth = .35, size = .7) +
    ggplot2::facet_wrap(~metric, nrow = 1) +
    ggplot2::scale_color_manual(values = palette_model, labels = c(fresh_group_batch = "Group+batch", fresh_MuSiC_PC1_5 = "MuSiC PC1-5")) +
    ggplot2::labs(x = NULL, y = "Point estimate and 95% percentile interval", color = "Discovery effects") +
    theme_publication() + ggplot2::theme(legend.position = "bottom", axis.text.x = ggplot2::element_text(size = 6, angle = 20, hjust = 1),
      strip.text = ggplot2::element_text(family = "Arial", face = "bold", size = 7))
  p16 <- p16a / p16b + patchwork::plot_layout(heights = c(1.15, .85)) +
    patchwork::plot_annotation(tag_levels = "a", title = "Fixed 346-gene effect agreement",
      subtitle = "Group-only cross-platform fits; 500 fixed stratified draws (seed 20261015); gate-failing estimable draws retained; shared four-control resampling",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(family = "Arial", face = "bold", size = 10),
        plot.subtitle = ggplot2::element_text(family = "Arial", size = 7),
        plot.tag = ggplot2::element_text(family = "Arial", face = "bold", size = 10)))
  write_plot_formats(p16, file.path(FIGURES, "Figure_S16_subgroup_agreement"), 11.5, 7.5)
  log_note("S15/S16 exported in PDF, PNG (600 dpi), and SVG")
  invisible(list(S15 = p15, S16 = p16))
}

if (identical(environment(), globalenv())) {
  args <- commandArgs(trailingOnly = TRUE)
  if ("--figures-only" %in% args) {
    log_note("Figure-only refresh requested; reusing locked fitted outputs")
    hubs <- as.character(read_csv(file.path(ANALYSIS, "locked_submitted_hub_set_346.csv"))$gene)
    baseline <- read_csv(file.path(ANALYSIS, "bulk_gene_level_baseline_full.csv.gz"))
    figure_state <- make_figures(list(hubs = hubs, fit_tables = list(baseline_full = baseline)), NULL)
  } else {
    bulk_state <- run_main()
    agreement_state <- run_agreement(bulk_state)
    figure_state <- make_figures(bulk_state, agreement_state)
  }
}

