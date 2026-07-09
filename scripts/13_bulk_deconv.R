#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

root <- Sys.getenv("PROJ_ROOT", getwd())
outdir <- file.path(root, "results", "deconv")
docdir <- file.path(root, "docs")

to_linear <- function(expr) {
  x <- as.matrix(expr)
  storage.mode(x) <- "double"
  shift <- min(x, na.rm = TRUE)
  y <- 2^(x - shift)
  y[!is.finite(y) | y <= 0] <- .Machine$double.eps
  dimnames(y) <- dimnames(x)
  y
}

run_nnls <- function(reference, bulk) {
  common <- intersect(rownames(reference), rownames(bulk))
  if (length(common) < 2L) stop("NNLS requires at least two common genes")
  ref <- as.matrix(reference[common, , drop = FALSE])
  bul <- as.matrix(bulk[common, , drop = FALSE])
  keep <- rowSums(is.finite(ref)) == ncol(ref) &
    apply(bul, 1, function(z) all(is.finite(z))) &
    rowSums(ref) > 0
  ref <- ref[keep, , drop = FALSE]
  bul <- bul[keep, , drop = FALSE]
  estimates <- t(vapply(seq_len(ncol(bul)), function(i) {
    fit <- nnls::nnls(ref, bul[, i])
    p <- pmax(stats::coef(fit), 0)
    if (!sum(p) > 0) p <- rep(1, length(p))
    p / sum(p)
  }, numeric(ncol(ref))))
  rownames(estimates) <- colnames(bul)
  colnames(estimates) <- colnames(ref)
  estimates
}

design_diagnostics <- function(meta, comp_scores = NULL) {
  dat <- meta
  if (!is.null(comp_scores) && ncol(comp_scores)) {
    dat <- cbind(dat, comp_scores)
  }
  rhs <- c("group", "batch", if (!is.null(comp_scores)) names(comp_scores) else NULL)
  design <- stats::model.matrix(stats::as.formula(paste("~", paste(rhs, collapse = "+"))),
                                data = dat)
  full_rank <- qr(design)$rank == ncol(design)
  group_col <- grep("^group", colnames(design), value = TRUE)[1]
  group_vif <- Inf
  if (!is.na(group_col) && full_rank) {
    others <- setdiff(colnames(design), c("(Intercept)", group_col))
    if (!length(others)) {
      group_vif <- 1
    } else {
      y <- design[, group_col]
      x <- design[, others, drop = FALSE]
      r2 <- summary(stats::lm(y ~ x))$r.squared
      group_vif <- if (is.finite(r2) && r2 < 1) 1 / (1 - r2) else Inf
    }
  }
  list(
    design = design,
    full_rank = full_rank,
    rank = qr(design)$rank,
    n_columns = ncol(design),
    condition_number = tryCatch(kappa(design), error = function(e) Inf),
    group_vif = group_vif
  )
}

read_expr <- function(path) {
  x <- read.csv(path, check.names = FALSE)
  genes <- x[[1]]
  m <- as.matrix(x[, -1, drop = FALSE])
  storage.mode(m) <- "double"
  rownames(m) <- genes
  m
}

write_fraction_table <- function(frac, meta, path) {
  stopifnot(all(rownames(frac) %in% meta$sample))
  out <- data.frame(sample = rownames(frac), frac, check.names = FALSE)
  out <- merge(meta, out, by = "sample", all.y = TRUE, sort = FALSE)
  out <- out[match(rownames(frac), out$sample), , drop = FALSE]
  write.csv(out, path, row.names = FALSE)
  out
}

select_reference_markers <- function(mean_log, n_per_type = 200L) {
  linear <- expm1(pmax(as.matrix(mean_log), 0))
  max_other <- sapply(seq_len(ncol(linear)), function(j) {
    apply(linear[, -j, drop = FALSE], 1, max)
  })
  specificity <- log2((linear + 0.05) / (max_other + 0.05))
  markers <- unique(unlist(lapply(seq_len(ncol(linear)), function(j) {
    eligible <- which(linear[, j] > 0.05 & is.finite(specificity[, j]))
    head(rownames(linear)[eligible[order(specificity[eligible, j], decreasing = TRUE)]],
         n_per_type)
  })))
  markers
}

validate_fractions <- function(x, expected_samples, label) {
  if (!identical(sort(rownames(x)), sort(expected_samples))) {
    stop(label, ": sample coverage mismatch")
  }
  if (any(!is.finite(x)) || any(x < -1e-8)) stop(label, ": invalid fractions")
  if (any(abs(rowSums(x) - 1) > 1e-5)) stop(label, ": fractions do not sum to one")
}

run_music <- function(bulk_linear, sce) {
  common <- intersect(rownames(bulk_linear), rownames(sce))
  if (length(common) < 1000L) stop("Too few common genes for MuSiC: ", length(common))
  ans <- MuSiC::music_prop(
    bulk.mtx = bulk_linear[common, , drop = FALSE],
    sc.sce = sce[common, ],
    clusters = "cell_type",
    samples = "sample",
    verbose = TRUE,
    normalize = TRUE
  )
  ans$Est.prop.weighted
}

group_tests <- function(frac_table, cohort, method, cell_types) {
  do.call(rbind, lapply(cell_types, function(ct) {
    control <- frac_table[frac_table$group == "Control", ct]
    tumor <- frac_table[frac_table$group == "Tumor", ct]
    wt <- tryCatch(stats::wilcox.test(tumor, control, exact = FALSE),
                   error = function(e) NULL)
    data.frame(
      cohort = cohort,
      method = method,
      cell_type = ct,
      n_control = length(control),
      n_tumor = length(tumor),
      median_control = median(control),
      median_tumor = median(tumor),
      median_shift = median(tumor) - median(control),
      p_value = if (is.null(wt)) NA_real_ else wt$p.value
    )
  }))
}

method_agreement <- function(music, nnls, cohort) {
  common_ct <- intersect(colnames(music), colnames(nnls))
  do.call(rbind, lapply(common_ct, function(ct) {
    data.frame(
      cohort = cohort,
      cell_type = ct,
      n_samples = nrow(music),
      spearman_r = suppressWarnings(stats::cor(music[, ct], nnls[, ct],
                                               method = "spearman"))
    )
  }))
}

significant_rows <- function(x, threshold = 0.05) {
  x[!is.na(x$padj) & x$padj < threshold, , drop = FALSE]
}

composition_pcs <- function(frac, max_pc = 5L) {
  z <- log(pmax(frac, 1e-6))
  clr <- z - rowMeans(z)
  pc <- stats::prcomp(clr, center = TRUE, scale. = FALSE)
  scores <- as.data.frame(pc$x[, seq_len(min(max_pc, ncol(pc$x))), drop = FALSE])
  rownames(scores) <- rownames(frac)
  scores
}

choose_adjustment <- function(meta, pc_scores, vif_max = 10, condition_max = 1000) {
  rows <- list()
  accepted <- character()
  for (pc in names(pc_scores)) {
    candidate <- c(accepted, pc)
    d <- design_diagnostics(meta, pc_scores[, candidate, drop = FALSE])
    pass <- d$full_rank && is.finite(d$group_vif) && d$group_vif <= vif_max &&
      is.finite(d$condition_number) && d$condition_number <= condition_max
    rows[[pc]] <- data.frame(
      candidate_pc = pc,
      included_pcs = paste(candidate, collapse = ";"),
      full_rank = d$full_rank,
      rank = d$rank,
      n_columns = d$n_columns,
      condition_number = d$condition_number,
      group_vif = d$group_vif,
      accepted = pass
    )
    if (pass) accepted <- candidate
  }
  list(accepted = accepted, diagnostics = do.call(rbind, rows))
}

fit_adjusted_deg <- function(expr, meta, pc_scores, accepted_pcs, original_deg) {
  dat <- meta
  dat$group <- factor(dat$group, levels = c("Control", "Tumor"))
  dat$batch <- factor(dat$batch)
  if (length(accepted_pcs)) dat <- cbind(dat, pc_scores[, accepted_pcs, drop = FALSE])
  rhs <- c("group", "batch", accepted_pcs)
  design <- model.matrix(as.formula(paste("~", paste(rhs, collapse = "+"))), data = dat)
  fit <- limma::eBayes(limma::lmFit(expr[, dat$sample, drop = FALSE], design),
                       trend = TRUE, robust = TRUE)
  tt <- limma::topTable(fit, coef = "groupTumor", number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  tt <- tt[, c("gene", setdiff(names(tt), "gene"))]
  names(tt)[names(tt) == "logFC"] <- "adjusted_logFC"
  names(tt)[names(tt) == "adj.P.Val"] <- "adjusted_adj.P.Val"
  names(tt)[names(tt) == "P.Value"] <- "adjusted_P.Value"
  old <- original_deg[, c("gene", "logFC", "adj.P.Val")]
  names(old)[-1] <- c("original_logFC", "original_adj.P.Val")
  out <- merge(old, tt, by = "gene", all = TRUE)
  out$original_sig <- !is.na(out$original_adj.P.Val) &
    out$original_adj.P.Val < 0.05 & abs(out$original_logFC) > 1
  out$adjusted_sig <- !is.na(out$adjusted_adj.P.Val) &
    out$adjusted_adj.P.Val < 0.05 & abs(out$adjusted_logFC) > 1
  out$classification <- ifelse(out$original_sig & out$adjusted_sig, "retained",
                               ifelse(out$original_sig & !out$adjusted_sig, "attenuated",
                                      ifelse(!out$original_sig & out$adjusted_sig,
                                             "adjusted-only", "neither")))
  out
}

plot_fractions <- function(tab, cell_types, title, path) {
  long <- do.call(rbind, lapply(cell_types, function(ct) {
    data.frame(sample = tab$sample, group = tab$group, cell_type = ct,
               fraction = tab[[ct]])
  }))
  p <- ggplot2::ggplot(long, ggplot2::aes(group, fraction, fill = group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.65) +
    ggplot2::geom_jitter(width = 0.12, size = 1.2, alpha = 0.75) +
    ggplot2::facet_wrap(~cell_type, scales = "free_y", ncol = 3) +
    ggplot2::scale_fill_manual(values = c(Control = "#2166AC", Tumor = "#B2182B")) +
    ggplot2::labs(title = title, x = NULL, y = "Estimated relative fraction") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "none")
  ggplot2::ggsave(path, p, width = 10, height = 7, dpi = 300)
}

main <- function() {
  suppressPackageStartupMessages({
    library(SeuratObject)
    library(SingleCellExperiment)
    library(MuSiC)
    library(nnls)
    library(limma)
    library(ggplot2)
  })
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  message("Reading persisted Seurat object")
  obj <- readRDS(file.path(root, "data", "processed", "GSE230375_seurat.rds"))
  counts <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
  meta_sc <- obj@meta.data[, c("sample", "cell_type"), drop = FALSE]
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(meta_sc)
  )
  rm(obj, counts)
  gc()

  mean_ref <- read_expr(file.path(root, "data", "processed",
                                  "GSE230375_celltype_meanexpr.csv"))
  markers <- select_reference_markers(mean_ref)
  writeLines(markers, file.path(outdir, "nnls_reference_markers.txt"))
  nnls_ref <- expm1(pmax(mean_ref[markers, , drop = FALSE], 0))

  cohorts <- list(
    GSE39645 = list(
      expr = read_expr(file.path(root, "data", "processed", "GSE39645_expr_gene.csv")),
      meta = read.csv(file.path(root, "data", "processed", "GSE39645_meta.csv"))
    ),
    GSE108524 = list(
      expr = read_expr(file.path(root, "data", "processed", "GSE108524_expr_gene.csv")),
      meta = read.csv(file.path(root, "data", "processed", "GSE108524_meta.csv"))
    )
  )
  cohorts$GSE108524$meta$group <- ifelse(cohorts$GSE108524$meta$group == "Control",
                                         "Control", "Tumor")

  all_tests <- list()
  all_agreement <- list()
  fraction_results <- list()
  for (cohort in names(cohorts)) {
    expr <- cohorts[[cohort]]$expr
    meta <- cohorts[[cohort]]$meta
    linear <- to_linear(expr)
    message("Running MuSiC: ", cohort)
    music <- run_music(linear, sce)
    music <- music[, colnames(mean_ref), drop = FALSE]
    message("Running marker NNLS: ", cohort)
    nnls_est <- run_nnls(nnls_ref, linear)
    nnls_est <- nnls_est[, colnames(mean_ref), drop = FALSE]
    validate_fractions(music, meta$sample, paste(cohort, "MuSiC"))
    validate_fractions(nnls_est, meta$sample, paste(cohort, "NNLS"))
    mt <- write_fraction_table(music, meta,
      file.path(outdir, paste0("cell_fractions_", cohort, "_music.csv")))
    nt <- write_fraction_table(nnls_est, meta,
      file.path(outdir, paste0("cell_fractions_", cohort, "_nnls.csv")))
    fraction_results[[cohort]] <- list(music = music, nnls = nnls_est,
                                       music_table = mt, nnls_table = nt)
    all_tests[[paste0(cohort, "_music")]] <- group_tests(
      mt, cohort, "MuSiC", colnames(mean_ref))
    all_tests[[paste0(cohort, "_nnls")]] <- group_tests(
      nt, cohort, "NNLS", colnames(mean_ref))
    all_agreement[[cohort]] <- method_agreement(music, nnls_est, cohort)
    plot_fractions(mt, colnames(mean_ref),
                   paste0(cohort, ": MuSiC relative fractions"),
                   file.path(outdir, paste0("fig_fractions_", cohort, ".png")))
  }
  tests <- do.call(rbind, all_tests)
  tests$padj <- ave(tests$p_value, interaction(tests$cohort, tests$method),
                    FUN = function(x) p.adjust(x, "BH"))
  write.csv(tests, file.path(outdir, "cell_fraction_group_tests.csv"), row.names = FALSE)
  agreement <- do.call(rbind, all_agreement)
  write.csv(agreement, file.path(outdir, "method_agreement.csv"), row.names = FALSE)

  expr396 <- cohorts$GSE39645$expr
  meta396 <- cohorts$GSE39645$meta
  meta396$group <- factor(meta396$group, levels = c("Control", "Tumor"))
  meta396$batch <- factor(meta396$batch)
  music396 <- fraction_results$GSE39645$music[meta396$sample, , drop = FALSE]
  pcs <- composition_pcs(music396)
  selection <- choose_adjustment(meta396, pcs)
  write.csv(selection$diagnostics, file.path(outdir, "design_diagnostics.csv"),
            row.names = FALSE)

  original <- read.csv(file.path(root, "results", "deg", "deg_full_gene.csv"))
  if (length(selection$accepted)) {
    adjusted <- fit_adjusted_deg(expr396, meta396, pcs, selection$accepted, original)
    write.csv(adjusted, file.path(outdir, "deg_composition_adjusted.csv"),
              row.names = FALSE)
    summary_tab <- as.data.frame(table(classification = adjusted$classification))
    write.csv(summary_tab, file.path(outdir, "deg_adjustment_summary.csv"),
              row.names = FALSE)
    plot_dat <- adjusted[adjusted$original_sig | adjusted$adjusted_sig, ]
    p <- ggplot2::ggplot(plot_dat,
                         ggplot2::aes(original_logFC, adjusted_logFC,
                                     color = classification)) +
      ggplot2::geom_hline(yintercept = 0, color = "grey70") +
      ggplot2::geom_vline(xintercept = 0, color = "grey70") +
      ggplot2::geom_point(alpha = 0.55, size = 1) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2) +
      ggplot2::labs(title = "GSE39645 DEG before and after composition adjustment",
                    subtitle = paste("Retained composition PCs:",
                                     paste(selection$accepted, collapse = ", ")),
                    x = "Original limma log2FC", y = "Composition-adjusted log2FC") +
      ggplot2::theme_bw(base_size = 11)
    ggplot2::ggsave(file.path(outdir, "fig_deg_adjustment.png"),
                    p, width = 7, height = 6, dpi = 300)
  } else {
    write.csv(data.frame(
      status = "not_estimable",
      reason = "No composition PC passed full-rank, condition-number and group-VIF gates"
    ), file.path(outdir, "deg_adjustment_summary.csv"), row.names = FALSE)
  }

  sig_tests <- significant_rows(tests)
  adjustment_text <- if (exists("adjusted")) {
    original_rows <- adjusted[adjusted$original_sig, , drop = FALSE]
    retained_n <- sum(original_rows$classification == "retained")
    attenuated_n <- sum(original_rows$classification == "attenuated")
    fc_r <- cor(original_rows$original_logFC, original_rows$adjusted_logFC,
                use = "complete.obs")
    direction_rate <- mean(sign(original_rows$original_logFC) ==
                             sign(original_rows$adjusted_logFC), na.rm = TRUE)
    c(
      paste0("- Original DEGs retained: ", retained_n, "/", nrow(original_rows),
             " (", round(100 * retained_n / nrow(original_rows), 1), "%)."),
      paste0("- Original DEGs attenuated: ", attenuated_n, "/", nrow(original_rows),
             " (", round(100 * attenuated_n / nrow(original_rows), 1), "%)."),
      paste0("- Original versus adjusted logFC: Pearson r=", round(fc_r, 3),
             "; direction agreement=", round(100 * direction_rate, 1), "%."),
      paste0("- Adjusted-only genes: ",
             sum(adjusted$classification == "adjusted-only"),
             "; treated as model-dependent reallocation, not new discoveries.")
    )
  } else {
    "- No composition-adjusted DEG model passed the estimability gate."
  }
  report <- c(
    "# 加深2: scRNA-reference bulk deconvolution and composition-adjusted DEG",
    "",
    "## Scope",
    "GSE230375 的六类粗粒度细胞作为参考，对 GSE39645 和 GSE108524 的 log2 微阵列表达进行 MuSiC 反卷积，并用 marker-restricted NNLS 作敏感性分析。微阵列先转换到正的线性尺度；估计值仅解释为跨平台 in-silico 相对比例分数，不等同于组织学细胞计数。",
    "",
    "## Key composition-adjustment findings",
    adjustment_text,
    "",
    "## Method agreement",
    paste(capture.output(print(agreement, row.names = FALSE)), collapse = "\n"),
    "",
    "## Group contrasts with FDR<0.05",
    if (nrow(sig_tests)) paste(capture.output(print(sig_tests, row.names = FALSE)),
                              collapse = "\n") else "No method/cohort contrast reached FDR<0.05.",
    "",
    "## Composition-adjustment gate",
    paste(capture.output(print(selection$diagnostics, row.names = FALSE)),
          collapse = "\n"),
    paste0("Accepted PCs: ",
           if (length(selection$accepted)) paste(selection$accepted, collapse = ", ")
           else "none; adjusted DEG was not considered estimable."),
    "",
    "## Interpretation constraints",
    "MuSiC 主要面向 bulk RNA-seq；本分析用于微阵列时属于跨平台探索性应用，并以 NNLS 方向一致性作敏感性检查。若组成轴与 VS 分组高度共线，则不能把调整后的 group 系数解释为独立的肿瘤内在效应。衰减仅表述为与组成依赖一致，不构成中介或因果证明。",
    "",
    "## Source paths",
    paste0("- Script: `", file.path(root, "scripts", "13_bulk_deconv.R"), "`"),
    paste0("- MuSiC/NNLS fractions and diagnostics: `", outdir, "/`"),
    paste0("- Single-cell reference: `",
           file.path(root, "data", "processed", "GSE230375_seurat.rds"), "`"),
    paste0("- Discovery expression/meta: `",
           file.path(root, "data", "processed", "GSE39645_expr_gene.csv"), "`, `",
           file.path(root, "data", "processed", "GSE39645_meta.csv"), "`"),
    paste0("- Validation expression/meta: `",
           file.path(root, "data", "processed", "GSE108524_expr_gene.csv"), "`, `",
           file.path(root, "data", "processed", "GSE108524_meta.csv"), "`")
  )
  writeLines(report, file.path(docdir, "deconv_report.md"))
  message("加深2完成")
}

if (sys.nframe() == 0L) main()
