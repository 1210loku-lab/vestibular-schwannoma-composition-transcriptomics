#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

ROOT <- Sys.getenv("VS_INPUT_ROOT", unset = ".")
EXPORT_ROOT <- Sys.getenv("VS_OUTPUT_ROOT", unset = ".")
ANALYSIS <- file.path(EXPORT_ROOT, "outputs")
SCRIPT_PATH <- file.path(EXPORT_ROOT, "scripts", "run_strong_simulation.R")
cell_types <- c("Schwann/肿瘤细胞", "巨噬/髓系", "T细胞", "成纤维", "内皮", "周细胞")
routes <- c("naive", "oracle_true_composition", "estimated_music_composition")
gene_classes <- c("null", "composition_only", "expression_only", "both", "reserved_marker", "non_target_background")

read_csv <- function(path) utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
read_gz_csv <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  utils::read.csv(con, stringsAsFactors = FALSE, check.names = FALSE)
}
write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "", quote = TRUE)
write_gz_csv <- function(x, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.csv(x, con, row.names = FALSE, na = "", quote = TRUE)
}
sha256 <- function(path) {
  out <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  if (length(out) != 1L || !grepl("^[0-9a-fA-F]{64}\\s", out)) stop("SHA-256 failed: ", path)
  strsplit(out, "\\s+")[[1]][1]
}

manifest_path <- file.path(ANALYSIS, "simulation_strong_FROZEN_INPUT_MANIFEST.csv")
frozen_manifest <- read_csv(manifest_path)
if (!all(file.exists(frozen_manifest$path))) stop("A frozen simulation input is missing")
rehashed <- vapply(frozen_manifest$path, sha256, character(1))
if (any(rehashed != frozen_manifest$sha256)) stop("Frozen simulation input hash mismatch")

required <- c("MuSiC", "SingleCellExperiment", "S4Vectors", "limma", "data.table")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Required package(s) missing: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(S4Vectors)
})

message("simulation strong extension execution started; protocol, generated inputs, and routes were frozen first.")
gene_truth <- read_csv(file.path(ANALYSIS, "simulation_gene_truth_and_marker_allocation.csv"))
if (nrow(gene_truth) != 1500L || anyDuplicated(gene_truth$gene)) stop("Invalid frozen gene truth table")
if (sum(gene_truth$reserved_marker) != 300L || any(gene_truth$reserved_marker & gene_truth$gene_class != "reserved_marker")) {
  stop("Reserved-marker allocation failed integrity check")
}
markers <- gene_truth$gene[gene_truth$reserved_marker]
if (length(intersect(markers, gene_truth$gene[gene_truth$gene_class %in% c("null", "composition_only", "expression_only", "both")]))) {
  stop("Reference marker identifiers overlap target evaluation genes")
}

ref_tab <- read_gz_csv(file.path(ANALYSIS, "simulation_synthetic_reference_counts.csv.gz"))
ref_genes <- as.character(ref_tab[[1]])
ref_counts <- as.matrix(ref_tab[, -1, drop = FALSE])
storage.mode(ref_counts) <- "double"
rownames(ref_counts) <- ref_genes
ref_meta <- read_csv(file.path(ANALYSIS, "simulation_synthetic_reference_cell_metadata.csv"))
if (ncol(ref_counts) != 900L || nrow(ref_meta) != ncol(ref_counts)) stop("Synthetic reference cell count mismatch")
if (!identical(colnames(ref_counts), ref_meta$cell)) stop("Reference count columns do not align with cell metadata")
if (!setequal(unique(ref_meta$cell_type), cell_types) || length(unique(ref_meta$sample)) != 5L) {
  stop("Synthetic reference donor or compartment coverage mismatch")
}
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = ref_counts),
  colData = S4Vectors::DataFrame(ref_meta[, c("sample", "cell_type"), drop = FALSE])
)

mix_tab <- read_gz_csv(file.path(ANALYSIS, "simulation_strong_mixture_log2_expression.csv.gz"))
mixture_genes <- as.character(mix_tab[[1]])
mixture_y <- as.matrix(mix_tab[, -1, drop = FALSE])
storage.mode(mixture_y) <- "double"
rownames(mixture_y) <- mixture_genes
sample_meta <- read_csv(file.path(ANALYSIS, "simulation_strong_sample_metadata.csv"))
composition_truth <- read_gz_csv(file.path(ANALYSIS, "simulation_strong_fraction_truth.csv.gz"))
if (nrow(sample_meta) != 800L || ncol(mixture_y) != nrow(sample_meta)) stop("Frozen mixture sample count mismatch")
if (!identical(colnames(mixture_y), sample_meta$sample)) stop("Mixture columns do not align with sample metadata")
if (anyDuplicated(sample_meta$sample) || anyDuplicated(sample_meta[c("scenario", "replicate", "sample")])) stop("Duplicate simulated sample keys")
if (nrow(composition_truth) != 4800L) stop("Frozen cell/RNA truth table has an unexpected row count")

fraction_matrix <- function(sample_ids, fraction_column) {
  out <- matrix(NA_real_, nrow = length(sample_ids), ncol = length(cell_types),
                dimnames = list(sample_ids, cell_types))
  for (i in seq_along(sample_ids)) {
    rows <- composition_truth[composition_truth$sample == sample_ids[i], , drop = FALSE]
    out[i, ] <- rows[[fraction_column]][match(cell_types, rows$cell_type)]
  }
  if (any(!is.finite(out))) stop("Truth fraction matrix is incomplete: ", fraction_column)
  out
}

design_diagnostics <- function(meta, pc_scores = NULL) {
  dat <- data.frame(group = factor(meta$group, levels = c("Control", "Tumor")))
  if (!is.null(pc_scores)) dat <- cbind(dat, as.data.frame(pc_scores))
  design <- stats::model.matrix(~ ., data = dat)
  rank <- qr(design)$rank
  group_col <- match("groupTumor", colnames(design))
  group_vif <- NA_real_
  if (!is.na(group_col) && rank == ncol(design)) {
    others <- design[, -group_col, drop = FALSE]
    y <- design[, group_col]
    fit_vif <- stats::lm.fit(others, y)
    tss <- sum((y - mean(y))^2)
    r2 <- if (tss > 0) 1 - sum(fit_vif$residuals^2) / tss else 1
    group_vif <- if (is.finite(r2) && r2 < 1) 1 / (1 - r2) else Inf
  }
  cond <- tryCatch(kappa(design), error = function(e) Inf)
  list(
    design = design,
    full_rank = rank == ncol(design),
    rank = rank,
    n_columns = ncol(design),
    condition_number = cond,
    group_vif = group_vif,
    gate_pass = rank == ncol(design) && is.finite(cond) && cond <= 1000 &&
      is.finite(group_vif) && group_vif <= 10
  )
}

composition_pcs <- function(frac) {
  z <- log(pmax(frac, 1e-6))
  clr <- z - rowMeans(z)
  pc <- stats::prcomp(clr, center = TRUE, scale. = FALSE)
  if (ncol(pc$x) < 5L) stop("Fewer than five composition PCs are available")
  scores <- pc$x[, seq_len(5L), drop = FALSE]
  colnames(scores) <- paste0("PC", seq_len(5L))
  scores
}

music_log_path <- file.path(ANALYSIS, "simulation_strong_MUSIC_RUN.log")
music_log <- file(music_log_path, open = "wt")
writeLines(c(
  paste0("simulation strong MuSiC log; package version ", as.character(utils::packageVersion("MuSiC"))),
  "Call: MuSiC::music_prop with explicit reserved unaffected markers, clusters=cell_type, samples=sample, normalize=TRUE, iter.max=1000, nu=1e-4, eps=0.01, centered=FALSE."
), music_log)

fit_one_route <- function(expr, meta, pc_scores, route, scenario_name, replicate_i, warning_rows) {
  scores <- if (route == "naive") NULL else pc_scores
  diag <- design_diagnostics(meta, scores)
  base_gene <- data.frame(
    scenario = scenario_name,
    replicate = replicate_i,
    route = route,
    gene = gene_truth$gene,
    gene_class = gene_truth$gene_class,
    beta = gene_truth$beta,
    logFC = NA_real_,
    raw_se = NA_real_,
    moderated_se = NA_real_,
    p_value = NA_real_,
    fdr = NA_real_,
    call_fdr05_abs_gt1 = NA,
    direction_error = NA,
    coefficient_error = NA_real_,
    fit_status = if (diag$full_rank) "not_run" else "rank_deficient",
    stringsAsFactors = FALSE
  )
  if (!diag$full_rank) return(list(genes = base_gene, diag = diag, warnings = warning_rows, fit_status = "rank_deficient_not_fitted"))
  dat <- data.frame(group = factor(meta$group, levels = c("Control", "Tumor")))
  if (!is.null(scores)) dat <- cbind(dat, as.data.frame(scores))
  design <- diag$design
  warning_text <- character()
  outcome <- tryCatch({
    lmfit <- withCallingHandlers(limma::lmFit(expr, design), warning = function(w) {
      warning_text <<- c(warning_text, paste0("lmFit: ", conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
    eb <- withCallingHandlers(limma::eBayes(lmfit, trend = TRUE, robust = TRUE), warning = function(w) {
      warning_text <<- c(warning_text, paste0("eBayes: ", conditionMessage(w)))
      invokeRestart("muffleWarning")
    })
    list(lmfit = lmfit, eb = eb)
  }, error = function(e) e)
  for (warning_msg in warning_text) {
    warning_rows[[length(warning_rows) + 1L]] <- data.frame(
      scenario = scenario_name, replicate = replicate_i, route = route,
      warning_type = "limma_warning", message = warning_msg, stringsAsFactors = FALSE
    )
  }
  if (inherits(outcome, "error")) {
    base_gene$fit_status <- paste0("fit_error: ", conditionMessage(outcome))
    return(list(genes = base_gene, diag = diag, warnings = warning_rows, fit_status = "fit_error"))
  }
  lmfit <- outcome$lmfit
  eb <- outcome$eb
  coef_name <- "groupTumor"
  if (!(coef_name %in% colnames(eb$coefficients))) {
    base_gene$fit_status <- "group_coefficient_missing"
    return(list(genes = base_gene, diag = diag, warnings = warning_rows, fit_status = "group_coefficient_missing"))
  }
  tt <- limma::topTable(eb, coef = coef_name, number = Inf, sort.by = "none")
  row_index <- match(gene_truth$gene, rownames(tt))
  coef_index <- match(coef_name, colnames(lmfit$coefficients))
  se_index <- match(coef_name, colnames(lmfit$stdev.unscaled))
  logfc <- as.numeric(eb$coefficients[row_index, coef_index])
  raw_se <- as.numeric(lmfit$stdev.unscaled[row_index, se_index] * lmfit$sigma[row_index])
  mod_se <- as.numeric(eb$stdev.unscaled[row_index, se_index] * sqrt(eb$s2.post[row_index]))
  pval <- as.numeric(tt$P.Value[row_index])
  fdr <- as.numeric(tt$adj.P.Val[row_index])
  call <- !is.na(fdr) & fdr < 0.05 & abs(logfc) > 1
  base_gene$logFC <- logfc
  base_gene$raw_se <- raw_se
  base_gene$moderated_se <- mod_se
  base_gene$p_value <- pval
  base_gene$fdr <- fdr
  base_gene$call_fdr05_abs_gt1 <- call
  base_gene$direction_error <- gene_truth$beta != 0 & sign(logfc) != sign(gene_truth$beta)
  base_gene$coefficient_error <- logfc - gene_truth$beta
  base_gene$fit_status <- "fit_ok"
  list(genes = base_gene, diag = diag, warnings = warning_rows, fit_status = "fit_ok")
}

fit_rows <- list()
diag_rows <- list()
warning_rows <- list()
repeat_rows <- list()
resource_rows <- list()
fraction_rows <- list()
fraction_metric_rows <- list()
repeat_count <- 0L
all_keys <- unique(sample_meta[, c("scenario", "replicate")])
total_timer <- proc.time()[[3]]

for (key_i in seq_len(nrow(all_keys))) {
  scenario_name <- all_keys$scenario[key_i]
  replicate_i <- all_keys$replicate[key_i]
  row_ids <- which(sample_meta$scenario == scenario_name & sample_meta$replicate == replicate_i)
  meta <- sample_meta[row_ids, , drop = FALSE]
  sample_ids <- meta$sample
  expr <- mixture_y[, sample_ids, drop = FALSE]
  frac_true <- fraction_matrix(sample_ids, "cell_fraction")
  frac_rna <- fraction_matrix(sample_ids, "rna_contribution_fraction")
  start <- proc.time()[[3]]
  gc(reset = TRUE)

  music_status <- "ok"
  music_error <- ""
  music_messages <- character()
  music_warnings <- character()
  music_fit <- tryCatch({
    withCallingHandlers(
      MuSiC::music_prop(
        bulk.mtx = 2^expr,
        sc.sce = sce,
        markers = markers,
        clusters = "cell_type",
        samples = "sample",
        verbose = TRUE,
        iter.max = 1000,
        nu = 1e-4,
        eps = 0.01,
        centered = FALSE,
        normalize = TRUE
      ),
      message = function(m) {
        music_messages <<- c(music_messages, conditionMessage(m))
        invokeRestart("muffleMessage")
      },
      warning = function(w) {
        music_warnings <<- c(music_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  }, error = function(e) e)
  if (inherits(music_fit, "error")) {
    music_status <- "error"
    music_error <- conditionMessage(music_fit)
    frac_est <- NULL
  } else {
    frac_est <- tryCatch(as.matrix(music_fit$Est.prop.weighted), error = function(e) NULL)
    if (is.null(frac_est) || nrow(frac_est) != length(sample_ids) ||
        !all(cell_types %in% colnames(frac_est)) || is.null(rownames(frac_est)) ||
        !setequal(rownames(frac_est), sample_ids)) {
      music_status <- "invalid_output_shape_or_labels"
      music_error <- "Est.prop.weighted does not match the frozen sample/compartment keys"
      frac_est <- NULL
    } else {
      frac_est <- frac_est[sample_ids, cell_types, drop = FALSE]
      if (any(!is.finite(frac_est)) || any(frac_est < -1e-8) || any(rowSums(frac_est) <= 0)) {
        music_status <- "invalid_fraction_values"
        music_error <- "Est.prop.weighted has non-finite, negative, or zero-sum rows"
        frac_est <- NULL
      } else {
        frac_est[frac_est < 0] <- 0
        frac_est <- frac_est / rowSums(frac_est)
      }
    }
  }
  log_lines <- c(
    paste0("\n[", scenario_name, " replicate=", replicate_i, "] status=", music_status),
    paste0("message: ", music_messages),
    paste0("warning: ", music_warnings),
    if (nzchar(music_error)) paste0("error: ", music_error) else character()
  )
  writeLines(log_lines, music_log)
  if (length(music_warnings)) {
    for (warning_msg in music_warnings) warning_rows[[length(warning_rows) + 1L]] <- data.frame(
      scenario = scenario_name, replicate = replicate_i, route = "music_prop",
      warning_type = "MuSiC_warning", message = warning_msg, stringsAsFactors = FALSE
    )
  }
  if (music_status != "ok") warning_rows[[length(warning_rows) + 1L]] <- data.frame(
    scenario = scenario_name, replicate = replicate_i, route = "music_prop",
    warning_type = paste0("MuSiC_", music_status), message = music_error, stringsAsFactors = FALSE
  )
  if (!is.null(frac_est)) {
    for (sample_i in seq_along(sample_ids)) {
      sample_id <- sample_ids[sample_i]
      for (ct_i in seq_along(cell_types)) {
        fraction_rows[[length(fraction_rows) + 1L]] <- data.frame(
          scenario = scenario_name,
          replicate = replicate_i,
          sample = sample_id,
          group = meta$group[sample_i],
          cell_type = cell_types[ct_i],
          cell_fraction_truth = frac_true[sample_i, ct_i],
          rna_contribution_truth = frac_rna[sample_i, ct_i],
          music_estimated_fraction = frac_est[sample_i, ct_i],
          abs_error_vs_cell_fraction = abs(frac_est[sample_i, ct_i] - frac_true[sample_i, ct_i]),
          abs_error_vs_rna_contribution = abs(frac_est[sample_i, ct_i] - frac_rna[sample_i, ct_i]),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  pcs_true <- tryCatch(composition_pcs(frac_true), error = function(e) e)
  pcs_music <- if (is.null(frac_est)) simpleError(music_error) else tryCatch(composition_pcs(frac_est), error = function(e) e)
  route_pcs <- list(naive = NULL, oracle_true_composition = pcs_true, estimated_music_composition = pcs_music)
  route_result <- list()
  for (route in routes) {
    pc_obj <- route_pcs[[route]]
    if (route != "naive" && inherits(pc_obj, "error")) {
      diag_stub <- list(full_rank = NA, rank = NA_integer_, n_columns = NA_integer_,
                        condition_number = NA_real_, group_vif = NA_real_, gate_pass = NA)
      gene_stub <- data.frame(
        scenario = scenario_name, replicate = replicate_i, route = route,
        gene = gene_truth$gene, gene_class = gene_truth$gene_class, beta = gene_truth$beta,
        logFC = NA_real_, raw_se = NA_real_, moderated_se = NA_real_, p_value = NA_real_,
        fdr = NA_real_, call_fdr05_abs_gt1 = NA, direction_error = NA,
        coefficient_error = NA_real_, fit_status = paste0("composition_PC_error: ", conditionMessage(pc_obj)),
        stringsAsFactors = FALSE
      )
      route_result[[route]] <- list(genes = gene_stub, diag = diag_stub,
                                    warnings = warning_rows, fit_status = "composition_PC_error")
    } else {
      route_result[[route]] <- fit_one_route(expr, meta, pc_obj, route, scenario_name,
                                             replicate_i, warning_rows)
      warning_rows <- route_result[[route]]$warnings
    }
    d <- route_result[[route]]$diag
    diag_rows[[length(diag_rows) + 1L]] <- data.frame(
      scenario = scenario_name,
      replicate = replicate_i,
      route = route,
      n_samples = nrow(meta),
      full_rank = d$full_rank,
      rank = d$rank,
      n_columns = d$n_columns,
      condition_number_unscaled = d$condition_number,
      group_vif = d$group_vif,
      gate_full_rank = d$full_rank,
      gate_condition_le_1000 = is.finite(d$condition_number) && d$condition_number <= 1000,
      gate_vif_le_10 = is.finite(d$group_vif) && d$group_vif <= 10,
      all_original_gates_pass = d$gate_pass,
      fit_status = route_result[[route]]$fit_status,
      music_status = if (route == "estimated_music_composition") music_status else "not_applicable",
      stringsAsFactors = FALSE
    )
    fit_rows[[length(fit_rows) + 1L]] <- route_result[[route]]$genes

    genes <- route_result[[route]]$genes
    for (class_name in gene_classes) {
      subset <- genes[genes$gene_class == class_name, , drop = FALSE]
      effect <- subset$beta != 0
      null <- !effect
      valid_calls <- !is.na(subset$call_fdr05_abs_gt1)
      calls <- subset$call_fdr05_abs_gt1
      repeat_rows[[length(repeat_rows) + 1L]] <- data.frame(
        scenario = scenario_name,
        replicate = replicate_i,
        route = route,
        gene_class = class_name,
        n_genes = nrow(subset),
        n_nonzero_beta = sum(effect),
        n_beta_zero = sum(null),
        n_fit_rows = sum(is.finite(subset$logFC)),
        threshold_call_rate_all = if (any(valid_calls)) mean(calls[valid_calls]) else NA_real_,
        false_call_rate_beta_zero = if (any(null & valid_calls)) mean(calls[null & valid_calls]) else NA_real_,
        sensitivity_nonzero_beta = if (any(effect & valid_calls)) mean(calls[effect & valid_calls]) else NA_real_,
        direction_error_rate_nonzero_beta = if (any(effect & is.finite(subset$logFC))) {
          mean(subset$direction_error[effect & is.finite(subset$logFC)], na.rm = TRUE)
        } else NA_real_,
        mean_coefficient_bias = if (any(is.finite(subset$coefficient_error))) mean(subset$coefficient_error, na.rm = TRUE) else NA_real_,
        mean_absolute_coefficient_error = if (any(is.finite(subset$coefficient_error))) mean(abs(subset$coefficient_error), na.rm = TRUE) else NA_real_,
        coefficient_rmse = if (any(is.finite(subset$coefficient_error))) sqrt(mean(subset$coefficient_error^2, na.rm = TRUE)) else NA_real_,
        median_raw_se = if (any(is.finite(subset$raw_se))) median(subset$raw_se, na.rm = TRUE) else NA_real_,
        median_moderated_se = if (any(is.finite(subset$moderated_se))) median(subset$moderated_se, na.rm = TRUE) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!is.null(frac_est)) {
    frac_tab <- do.call(rbind, fraction_rows[(length(fraction_rows) - length(sample_ids) * length(cell_types) + 1L):length(fraction_rows)])
    fraction_metric_rows[[length(fraction_metric_rows) + 1L]] <- data.frame(
      scenario = scenario_name,
      replicate = replicate_i,
      n_samples = nrow(meta),
      mean_absolute_error_vs_cell_fraction = mean(frac_tab$abs_error_vs_cell_fraction),
      rmse_vs_cell_fraction = sqrt(mean((frac_tab$music_estimated_fraction - frac_tab$cell_fraction_truth)^2)),
      mean_absolute_error_vs_rna_contribution = mean(frac_tab$abs_error_vs_rna_contribution),
      stringsAsFactors = FALSE
    )
  } else {
    fraction_metric_rows[[length(fraction_metric_rows) + 1L]] <- data.frame(
      scenario = scenario_name, replicate = replicate_i, n_samples = nrow(meta),
      mean_absolute_error_vs_cell_fraction = NA_real_, rmse_vs_cell_fraction = NA_real_,
      mean_absolute_error_vs_rna_contribution = NA_real_, stringsAsFactors = FALSE
    )
  }

  elapsed <- proc.time()[[3]] - start
  gc_observed <- gc()
  gc_peak_mb <- sum(gc_observed[, 7])
  repeat_count <- repeat_count + 1L
  resource_rows[[length(resource_rows) + 1L]] <- data.frame(
    scenario = scenario_name,
    replicate = replicate_i,
    route = "__resource__",
    gene_class = "__all__",
    n_genes = NA_integer_, n_nonzero_beta = NA_integer_, n_beta_zero = NA_integer_, n_fit_rows = NA_integer_,
    threshold_call_rate_all = NA_real_, false_call_rate_beta_zero = NA_real_, sensitivity_nonzero_beta = NA_real_,
    direction_error_rate_nonzero_beta = NA_real_, mean_coefficient_bias = NA_real_,
    mean_absolute_coefficient_error = NA_real_, coefficient_rmse = NA_real_, median_raw_se = NA_real_,
    median_moderated_se = NA_real_, runtime_seconds = elapsed, gc_peak_heap_mb_observed = gc_peak_mb,
    music_status = music_status, music_error = music_error,
    stringsAsFactors = FALSE
  )
  if (repeat_count %% 2L == 0L || repeat_count == 1L) {
    message(sprintf("Completed %d/20 repeats (%s, replicate %d); MuSiC=%s; %.2f sec.",
                    repeat_count, scenario_name, replicate_i, music_status, elapsed))
  }
}

total_elapsed <- proc.time()[[3]] - total_timer
close(music_log)
gene_fits <- do.call(rbind, fit_rows)
design_diag <- do.call(rbind, diag_rows)
warning_table <- if (length(warning_rows)) do.call(rbind, warning_rows) else data.frame(
  scenario = character(), replicate = integer(), route = character(), warning_type = character(), message = character()
)
metric_table <- do.call(rbind, repeat_rows)
resource_table <- do.call(rbind, resource_rows)
fraction_table <- if (length(fraction_rows)) do.call(rbind, fraction_rows) else data.frame(
  scenario = character(), replicate = integer(), sample = character(), group = character(), cell_type = character(),
  cell_fraction_truth = numeric(), rna_contribution_truth = numeric(), music_estimated_fraction = numeric(),
  abs_error_vs_cell_fraction = numeric(), abs_error_vs_rna_contribution = numeric()
)
fraction_metrics <- do.call(rbind, fraction_metric_rows)

# Pairwise effect and uncertainty changes from the naive route are retained per gene.
fit_key <- paste(gene_fits$scenario, gene_fits$replicate, gene_fits$gene, sep = "\r")
naive <- gene_fits[gene_fits$route == "naive", , drop = FALSE]
naive_key <- paste(naive$scenario, naive$replicate, naive$gene, sep = "\r")
comparison_rows <- list()
for (route in routes[routes != "naive"]) {
  adjusted <- gene_fits[gene_fits$route == route, , drop = FALSE]
  adjusted_key <- paste(adjusted$scenario, adjusted$replicate, adjusted$gene, sep = "\r")
  m <- match(adjusted_key, naive_key)
  comparison_rows[[route]] <- data.frame(
    scenario = adjusted$scenario,
    replicate = adjusted$replicate,
    route = route,
    gene = adjusted$gene,
    gene_class = adjusted$gene_class,
    delta_logFC_vs_naive = adjusted$logFC - naive$logFC[m],
    raw_se_ratio_vs_naive = adjusted$raw_se / naive$raw_se[m],
    moderated_se_ratio_vs_naive = adjusted$moderated_se / naive$moderated_se[m],
    stringsAsFactors = FALSE
  )
}
route_changes <- do.call(rbind, comparison_rows)
route_change_summary <- do.call(rbind, lapply(split(route_changes, interaction(route_changes$scenario,
                                                                                  route_changes$replicate,
                                                                                  route_changes$route,
                                                                                  route_changes$gene_class,
                                                                                  drop = TRUE)), function(x) {
  data.frame(
    scenario = x$scenario[1], replicate = x$replicate[1], route = x$route[1], gene_class = x$gene_class[1],
    median_abs_delta_logFC_vs_naive = if (any(is.finite(x$delta_logFC_vs_naive))) median(abs(x$delta_logFC_vs_naive), na.rm = TRUE) else NA_real_,
    median_raw_se_ratio_vs_naive = if (any(is.finite(x$raw_se_ratio_vs_naive))) median(x$raw_se_ratio_vs_naive, na.rm = TRUE) else NA_real_,
    median_moderated_se_ratio_vs_naive = if (any(is.finite(x$moderated_se_ratio_vs_naive))) median(x$moderated_se_ratio_vs_naive, na.rm = TRUE) else NA_real_,
    stringsAsFactors = FALSE
  )
}))

summarize_rate <- function(x) {
  vals <- x[is.finite(x)]
  if (!length(vals)) return(c(mean = NA_real_, mcse = NA_real_, sd = NA_real_, n = 0))
  c(mean = mean(vals), mcse = if (length(vals) > 1) stats::sd(vals) / sqrt(length(vals)) else NA_real_,
    sd = if (length(vals) > 1) stats::sd(vals) else NA_real_, n = length(vals))
}
metric_summary <- do.call(rbind, lapply(split(metric_table, interaction(metric_table$scenario,
                                                                         metric_table$route,
                                                                         metric_table$gene_class,
                                                                         drop = TRUE)), function(x) {
  call <- summarize_rate(x$threshold_call_rate_all)
  fp <- summarize_rate(x$false_call_rate_beta_zero)
  sens <- summarize_rate(x$sensitivity_nonzero_beta)
  dir <- summarize_rate(x$direction_error_rate_nonzero_beta)
  data.frame(
    scenario = x$scenario[1], route = x$route[1], gene_class = x$gene_class[1],
    attempted_repeats = 20L,
    repeats_with_fit = sum(x$n_fit_rows == x$n_genes),
    call_rate_mean = call["mean"], call_rate_mcse = call["mcse"], call_rate_sd = call["sd"],
    beta_zero_false_call_rate_mean = fp["mean"], beta_zero_false_call_rate_mcse = fp["mcse"],
    nonzero_beta_sensitivity_mean = sens["mean"], nonzero_beta_sensitivity_mcse = sens["mcse"],
    direction_error_rate_mean = dir["mean"], direction_error_rate_mcse = dir["mcse"],
    mean_coefficient_bias_across_repeats = mean(x$mean_coefficient_bias, na.rm = TRUE),
    mean_absolute_coefficient_error_across_repeats = mean(x$mean_absolute_coefficient_error, na.rm = TRUE),
    mean_coefficient_rmse_across_repeats = mean(x$coefficient_rmse, na.rm = TRUE),
    median_raw_se_across_repeats = median(x$median_raw_se, na.rm = TRUE),
    median_moderated_se_across_repeats = median(x$median_moderated_se, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

write_gz_csv(gene_fits, file.path(ANALYSIS, "simulation_strong_gene_fit_outputs.csv.gz"))
write_csv(design_diag, file.path(ANALYSIS, "simulation_strong_design_and_gate_diagnostics.csv"))
write_csv(warning_table, file.path(ANALYSIS, "simulation_strong_all_warnings_and_failures.csv"))
write_csv(metric_table, file.path(ANALYSIS, "simulation_strong_repeat_gene_class_metrics.csv"))
write_csv(metric_summary, file.path(ANALYSIS, "simulation_strong_pilot_metric_summary.csv"))
write_gz_csv(fraction_table, file.path(ANALYSIS, "simulation_strong_MuSiC_fraction_estimates_and_truth.csv.gz"))
write_csv(fraction_metrics, file.path(ANALYSIS, "simulation_strong_fraction_recovery_metrics_by_repeat.csv"))
write_gz_csv(route_changes, file.path(ANALYSIS, "simulation_strong_gene_effect_se_changes_vs_naive.csv.gz"))
write_csv(route_change_summary, file.path(ANALYSIS, "simulation_strong_effect_se_change_summary_by_repeat.csv"))
write_csv(resource_table[, c("scenario", "replicate", "runtime_seconds", "gc_peak_heap_mb_observed", "music_status", "music_error")],
          file.path(ANALYSIS, "simulation_strong_repeat_resource_metrics.csv"))

output_paths <- c(
  SCRIPT_PATH,
  file.path(ANALYSIS, "simulation_strong_FROZEN_INPUT_MANIFEST.csv"),
  file.path(ANALYSIS, "simulation_strong_MUSIC_RUN.log"),
  file.path(ANALYSIS, "simulation_strong_gene_fit_outputs.csv.gz"),
  file.path(ANALYSIS, "simulation_strong_design_and_gate_diagnostics.csv"),
  file.path(ANALYSIS, "simulation_strong_all_warnings_and_failures.csv"),
  file.path(ANALYSIS, "simulation_strong_repeat_gene_class_metrics.csv"),
  file.path(ANALYSIS, "simulation_strong_pilot_metric_summary.csv"),
  file.path(ANALYSIS, "simulation_strong_MuSiC_fraction_estimates_and_truth.csv.gz"),
  file.path(ANALYSIS, "simulation_strong_fraction_recovery_metrics_by_repeat.csv"),
  file.path(ANALYSIS, "simulation_strong_gene_effect_se_changes_vs_naive.csv.gz"),
  file.path(ANALYSIS, "simulation_strong_effect_se_change_summary_by_repeat.csv"),
  file.path(ANALYSIS, "simulation_strong_repeat_resource_metrics.csv")
)
out_manifest <- data.frame(
  path = normalizePath(output_paths, winslash = "/", mustWork = TRUE),
  bytes = file.info(output_paths)$size,
  sha256 = vapply(output_paths, sha256, character(1)),
  stringsAsFactors = FALSE
)
write_csv(out_manifest, file.path(ANALYSIS, "simulation_strong_OUTPUT_MANIFEST.csv"))

cat("simulation strong extension complete.\n")
cat("attempted scenario-repeats=", repeat_count, " (expected 20); samples=", nrow(sample_meta),
    "; routes=", paste(routes, collapse = ","), "; MuSiC non-ok=", sum(design_diag$route == "estimated_music_composition" & design_diag$music_status != "ok"),
    "; total_elapsed_seconds=", round(total_elapsed, 2), "\n", sep = "")
cat("design rows=", nrow(design_diag), "; gene-fit rows=", nrow(gene_fits),
    "; warning/failure records=", nrow(warning_table), "; fraction rows=", nrow(fraction_table), "\n", sep = "")
cat("maximum observed R vector heap across repeats (MB)=", round(max(resource_table$gc_peak_heap_mb_observed), 2), "\n", sep = "")
