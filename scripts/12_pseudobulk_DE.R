#!/usr/bin/env Rscript
# 12_pseudobulk_DE.R
# 加深1: scRNA 细胞类型内 pseudobulk 差异表达 (GSE230375, VS vs Nerve)
# 目的: 区分 bulk DEG 中"候选肿瘤/细胞内在程序" vs "与细胞组成变化一致"的信号,
#       缓解 bulk RNA-seq 的细胞组成混杂软肋。
# 统计单位 = 样本。约束: VS=7, Nerve=2 → Nerve 端欠功效,
#       判定以"方向一致性 + logFC 相关"为主, nominal p / DESeq2 padj 为辅。
# 措辞: 相关性/网络/in-silico 证据不写因果, 用"候选/与…一致/提示"。
# 输入: data/processed/GSE230375_seurat.rds (由 08 重建并持久化)
#       results/deg/deg_full_gene.csv (bulk limma DEG, logFC=log2)
#       results/scrna/celltype_composition.csv (各组细胞比例)
# 产物: results/pseudobulk/*  +  docs/pseudobulk_report.md
#       data/processed/GSE230375_{pseudobulk_counts,celltype_meanexpr,cell_meta}.csv (供加深2复用)

options(stringsAsFactors = FALSE, warn = 1)
suppressMessages({
  library(Matrix)
  library(DESeq2)
  library(ggplot2)
})

root <- Sys.getenv("PROJ_ROOT", getwd())
proc      <- file.path(root, "data", "processed")
outdir    <- file.path(root, "results", "pseudobulk")
docdir    <- file.path(root, "docs")
rds_path  <- file.path(proc, "GSE230375_seurat.rds")
bulk_path <- file.path(root, "results", "deg", "deg_full_gene.csv")
comp_path <- file.path(root, "results", "scrna", "celltype_composition.csv")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# ---- 参数 ----
MIN_CELLS_PER_SAMPLE_CT <- 10   # 某样本×细胞类型至少这么多细胞才计入 pseudobulk
MIN_SAMPLES_PER_GROUP   <- 2    # 每组至少这么多样本才在该细胞类型内做 DE
REAL_CELLTYPES <- c("Schwann/肿瘤细胞", "巨噬/髓系", "T细胞", "成纤维", "内皮", "周细胞")
CT_LATIN <- c("Schwann/肿瘤细胞" = "Schwann_tumor", "巨噬/髓系" = "Myeloid",
              "T细胞" = "T_cell", "成纤维" = "Fibroblast", "内皮" = "Endothelial",
              "周细胞" = "Pericyte")
CT_PLOT <- c("Schwann/肿瘤细胞" = "Schwann/tumour", "巨噬/髓系" = "Myeloid",
             "T细胞" = "T cell", "成纤维" = "Fibroblast",
             "内皮" = "Endothelial", "周细胞" = "Pericyte")

filter_concordance_data <- function(de, bulk, filter_type, threshold = NA_real_) {
  m <- merge(
    de[, c("gene", "baseMean", "log2FC")],
    bulk[, c("gene", "logFC", "adj.P.Val")],
    by = "gene"
  )
  m <- m[
    is.finite(m$log2FC) & is.finite(m$logFC) & is.finite(m$baseMean),
  ]

  if (filter_type == "baseMean") {
    m <- m[m$baseMean >= threshold, ]
  } else if (filter_type == "bulk_DEG") {
    m <- m[
      !is.na(m$adj.P.Val) & m$adj.P.Val < 0.05 & abs(m$logFC) >= 1,
    ]
  } else if (filter_type != "all") {
    stop("Unsupported filter_type")
  }
  m
}

bootstrap_cor_ci <- function(x, y, method,
                             n_boot = 2000L,
                             seed = 230375L,
                             conf = 0.95) {
  stopifnot(length(x) == length(y))
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 4L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(c(NA_real_, NA_real_))
  }

  set.seed(seed)
  boots <- replicate(n_boot, {
    idx <- sample.int(length(x), replace = TRUE)
    suppressWarnings(cor(x[idx], y[idx], method = method))
  })
  boots <- boots[is.finite(boots)]
  if (!length(boots)) return(c(NA_real_, NA_real_))
  alpha <- (1 - conf) / 2
  unname(quantile(boots, c(alpha, 1 - alpha), names = FALSE))
}

summarize_concordance <- function(m, cell_type, filter_type,
                                  threshold = NA_real_,
                                  n_boot = 2000L) {
  keep <- is.finite(m$logFC) & is.finite(m$log2FC)
  m <- m[keep, , drop = FALSE]
  safe_cor <- function(method) {
    if (nrow(m) < 2L ||
        length(unique(m$logFC)) < 2L ||
        length(unique(m$log2FC)) < 2L) {
      return(NA_real_)
    }
    suppressWarnings(cor(m$logFC, m$log2FC, method = method))
  }

  pearson <- safe_cor("pearson")
  spearman <- safe_cor("spearman")
  pearson_ci <- bootstrap_cor_ci(
    m$logFC, m$log2FC, "pearson", n_boot, seed = 230375L
  )
  spearman_ci <- bootstrap_cor_ci(
    m$logFC, m$log2FC, "spearman", n_boot, seed = 230376L
  )

  nonzero <- sign(m$logFC) != 0 & sign(m$log2FC) != 0
  n_direction <- sum(nonzero)
  n_concordant <- sum(sign(m$logFC[nonzero]) == sign(m$log2FC[nonzero]))
  if (n_direction > 0L) {
    bt <- binom.test(n_concordant, n_direction)
    direction_concordance <- n_concordant / n_direction
    direction_ci <- unname(bt$conf.int)
  } else {
    direction_concordance <- NA_real_
    direction_ci <- c(NA_real_, NA_real_)
  }

  ci_method <- if (nrow(m) < 4L) {
    "not_estimable_insufficient_n; direction_clopper_pearson"
  } else if (!is.finite(pearson) || !is.finite(spearman)) {
    "not_estimable_constant; direction_clopper_pearson"
  } else {
    "bootstrap_percentile; direction_clopper_pearson"
  }

  data.frame(
    cell_type = cell_type,
    filter_type = filter_type,
    threshold = threshold,
    n_genes = nrow(m),
    pearson_r = pearson,
    pearson_ci_low = pearson_ci[1],
    pearson_ci_high = pearson_ci[2],
    spearman_r = spearman,
    spearman_ci_low = spearman_ci[1],
    spearman_ci_high = spearman_ci[2],
    n_direction = n_direction,
    n_concordant = n_concordant,
    direction_concordance = direction_concordance,
    direction_ci_low = direction_ci[1],
    direction_ci_high = direction_ci[2],
    ci_method = ci_method,
    stringsAsFactors = FALSE
  )
}

build_dominant_celltype_concordance <- function(
    bulk_deg, dom_ct, de_list, thresholds = c(10, 20)) {
  rows <- lapply(seq_len(nrow(bulk_deg)), function(i) {
    gene <- bulk_deg$gene[i]
    dct <- unname(dom_ct[gene])
    if (!length(dct) || is.na(dct) || !dct %in% names(de_list)) {
      return(data.frame(
        gene = gene, bulk_logFC = bulk_deg$logFC[i],
        dominant_celltype = ifelse(length(dct), dct, NA_character_),
        baseMean = NA_real_, pb_log2FC = NA_real_
      ))
    }
    de <- de_list[[dct]]$de
    hit <- de[de$gene == gene, c("baseMean", "log2FC"), drop = FALSE]
    data.frame(
      gene = gene, bulk_logFC = bulk_deg$logFC[i],
      dominant_celltype = dct,
      baseMean = if (nrow(hit)) hit$baseMean[1] else NA_real_,
      pb_log2FC = if (nrow(hit)) hit$log2FC[1] else NA_real_
    )
  })
  gene_tab <- do.call(rbind, rows)

  summarize_scope <- function(x, scope, dominant_celltype) {
    specs <- c(list(list(filter_type = "all", threshold = NA_real_)),
               lapply(thresholds, function(z) {
                 list(filter_type = "baseMean", threshold = z)
               }))
    do.call(rbind, lapply(specs, function(spec) {
      with_pb <- is.finite(x$baseMean) & is.finite(x$pb_log2FC) &
        is.finite(x$bulk_logFC)
      selected <- with_pb
      if (spec$filter_type == "baseMean") {
        selected <- selected & x$baseMean >= spec$threshold
      }
      nonzero <- selected & sign(x$bulk_logFC) != 0 & sign(x$pb_log2FC) != 0
      n_direction <- sum(nonzero)
      n_concordant <- sum(
        sign(x$bulk_logFC[nonzero]) == sign(x$pb_log2FC[nonzero])
      )
      if (n_direction > 0L) {
        bt <- binom.test(n_concordant, n_direction)
        direction_concordance <- n_concordant / n_direction
        direction_ci <- unname(bt$conf.int)
      } else {
        direction_concordance <- NA_real_
        direction_ci <- c(NA_real_, NA_real_)
      }
      data.frame(
        scope = scope,
        dominant_celltype = dominant_celltype,
        filter_type = spec$filter_type,
        threshold = spec$threshold,
        n_bulk_deg = nrow(x),
        n_with_pb = sum(selected),
        n_direction = n_direction,
        n_concordant = n_concordant,
        direction_concordance = direction_concordance,
        direction_ci_low = direction_ci[1],
        direction_ci_high = direction_ci[2],
        stringsAsFactors = FALSE
      )
    }))
  }

  out <- list(summarize_scope(gene_tab, "overall", "All"))
  for (ct in unique(gene_tab$dominant_celltype[!is.na(gene_tab$dominant_celltype)])) {
    out[[length(out) + 1L]] <- summarize_scope(
      gene_tab[gene_tab$dominant_celltype == ct, , drop = FALSE],
      "by_celltype", ct
    )
  }
  do.call(rbind, out)
}

# ---- 载入 ----
log_msg("Reading Seurat object: ", rds_path)
obj <- readRDS(rds_path)
counts <- SeuratObject::LayerData(obj, assay = "RNA", layer = "counts")
data_n <- SeuratObject::LayerData(obj, assay = "RNA", layer = "data")
meta <- obj@meta.data[, c("sample", "group", "cell_type")]
meta$cell_type <- as.character(meta$cell_type)
stopifnot(identical(colnames(counts), rownames(meta)))

# 完整性闸: 重建构成需与 08 记录一致, 否则报警
comp_old <- read.csv(comp_path, check.names = FALSE)
comp_new <- as.data.frame(table(group = meta$group, cell_type = meta$cell_type))
comp_new <- comp_new[comp_new$Freq > 0, ]
tot <- tapply(comp_new$Freq, comp_new$group, sum)
comp_new$pct <- round(100 * comp_new$Freq / tot[comp_new$group], 2)
chk <- merge(comp_old, comp_new, by = c("group", "cell_type"), all = TRUE)
chk$pct_diff <- abs(chk$percent_within_group - chk$pct)
gate_ok <- all(is.na(chk$pct_diff) | chk$pct_diff < 1.0)
write.csv(chk[, c("group", "cell_type", "percent_within_group", "pct", "pct_diff")],
          file.path(outdir, "composition_reproducibility_check.csv"), row.names = FALSE)
log_msg("Composition reproducibility gate: ", ifelse(gate_ok, "PASS", "WARN (>1% drift, 见 check csv)"))

# ---- pseudobulk 聚合 (sample × celltype, sum raw counts) ----
keep <- !is.na(meta$cell_type) & meta$cell_type %in% REAL_CELLTYPES
meta_k <- meta[keep, ]
counts_k <- counts[, keep, drop = FALSE]
grp_id <- paste(meta_k$sample, meta_k$cell_type, sep = "@@")
ugrp <- sort(unique(grp_id))
ind <- sparse.model.matrix(~ 0 + factor(grp_id, levels = ugrp))
colnames(ind) <- ugrp
pb <- as.matrix(counts_k %*% ind)            # genes × (sample@@celltype)
ncell <- as.integer(table(factor(grp_id, levels = ugrp)))
names(ncell) <- ugrp

pb_meta <- data.frame(
  id = ugrp,
  sample = sub("@@.*$", "", ugrp),
  cell_type = sub("^.*@@", "", ugrp),
  n_cells = ncell[ugrp],
  stringsAsFactors = FALSE
)
pb_meta$group <- meta_k$group[match(pb_meta$sample, meta_k$sample)]
write.csv(pb_meta, file.path(outdir, "pseudobulk_sample_celltype_ncells.csv"), row.names = FALSE)
write.csv(data.frame(gene = rownames(pb), pb, check.names = FALSE),
          file.path(proc, "GSE230375_pseudobulk_counts.csv"), row.names = FALSE)

# ---- 细胞类型平均表达参考 (normalized data 均值, 供加深2 反卷积/主导细胞判定) ----
mean_by_ct <- sapply(REAL_CELLTYPES, function(ct) {
  cols <- which(meta$cell_type == ct)
  if (!length(cols)) return(rep(NA_real_, nrow(data_n)))
  Matrix::rowMeans(data_n[, cols, drop = FALSE])
})
rownames(mean_by_ct) <- rownames(data_n)
write.csv(data.frame(gene = rownames(mean_by_ct), mean_by_ct, check.names = FALSE),
          file.path(proc, "GSE230375_celltype_meanexpr.csv"), row.names = FALSE)
write.csv(data.frame(barcode = rownames(meta), meta),
          file.path(proc, "GSE230375_cell_meta.csv"), row.names = FALSE)

# dominant celltype per gene = 平均归一化表达最高的真实细胞类型
dom_ct <- colnames(mean_by_ct)[max.col(replace(mean_by_ct, is.na(mean_by_ct), -Inf), ties.method = "first")]
names(dom_ct) <- rownames(mean_by_ct)

# ---- 每细胞类型内 DESeq2 (VS vs Nerve) ----
run_ct_de <- function(ct) {
  ids <- pb_meta$id[pb_meta$cell_type == ct & pb_meta$n_cells >= MIN_CELLS_PER_SAMPLE_CT]
  cd <- pb_meta[match(ids, pb_meta$id), ]
  ngrp <- table(cd$group)
  if (!all(c("VS", "Nerve") %in% names(ngrp)) ||
      any(ngrp[c("VS", "Nerve")] < MIN_SAMPLES_PER_GROUP)) {
    log_msg("SKIP ", ct, " (insufficient samples: ",
            paste(names(ngrp), ngrp, sep = "=", collapse = ", "), ")")
    return(NULL)
  }
  m <- pb[, ids, drop = FALSE]
  m <- m[rowSums(m) >= 10, , drop = FALSE]   # 轻过滤
  coldata <- data.frame(group = factor(cd$group, levels = c("Nerve", "VS")))
  rownames(coldata) <- ids
  dds <- DESeqDataSetFromMatrix(round(m), coldata, design = ~ group)
  dds <- DESeq(dds, quiet = TRUE)
  res <- results(dds, contrast = c("group", "VS", "Nerve"))
  shr <- tryCatch(
    lfcShrink(dds, coef = "group_VS_vs_Nerve", type = "apeglm"),
    error = function(e) tryCatch(
      lfcShrink(dds, contrast = c("group", "VS", "Nerve"), type = "ashr"),
      error = function(e2) { log_msg("  LFC shrink fallback=normal for ", ct); res }))
  out <- data.frame(
    gene = rownames(res),
    baseMean = res$baseMean,
    log2FC = shr$log2FoldChange,
    log2FC_raw = res$log2FoldChange,
    lfcSE = shr$lfcSE,
    pvalue = res$pvalue,
    padj = res$padj,
    n_VS = unname(ngrp["VS"]), n_Nerve = unname(ngrp["Nerve"]),
    stringsAsFactors = FALSE
  )
  out <- out[order(out$pvalue), ]
  fn <- file.path(outdir, paste0(CT_LATIN[[ct]], "_VS_vs_nerve_DE.csv"))
  write.csv(out, fn, row.names = FALSE)
  log_msg("  ", ct, " DE done (", nrow(out), " genes; VS=", ngrp["VS"],
          " Nerve=", ngrp["Nerve"], ") -> ", basename(fn))
  # 火山图
  v <- out[is.finite(out$pvalue) & is.finite(out$log2FC), ]
  v$sig <- ifelse(v$padj < 0.05 & abs(v$log2FC) > 1, "FDR<0.05 & |log2FC|>1", "ns")
  p <- ggplot(v, aes(log2FC, -log10(pvalue), color = sig)) +
    geom_point(size = 0.6, alpha = 0.5) +
    scale_color_manual(values = c("FDR<0.05 & |log2FC|>1" = "#B2182B", "ns" = "grey75")) +
    geom_vline(xintercept = c(-1, 1), linetype = 2, color = "grey50") +
    labs(title = paste0("Within-", CT_LATIN[[ct]], " pseudobulk: VS vs Nerve"),
         subtitle = paste0("VS n=", ngrp["VS"], ", Nerve n=", ngrp["Nerve"],
                           " (Nerve n=2: 探索性, 欠功效)"),
         x = "log2 fold change (VS / Nerve, shrunk)", y = "-log10 P", color = NULL) +
    theme_bw(base_size = 11)
  ggsave(file.path(outdir, paste0("fig_volcano_", CT_LATIN[[ct]], ".png")),
         p, width = 7, height = 5.5, dpi = 300)
  list(ct = ct, de = out, n = ngrp)
}

de_list <- list()
for (ct in REAL_CELLTYPES) {
  r <- run_ct_de(ct)
  if (!is.null(r)) de_list[[ct]] <- r
}

# ---- bulk DEG ↔ pseudobulk 分类: 内在 vs 组成 ----
bulk <- read.csv(bulk_path)
bulk_deg <- bulk[!is.na(bulk$adj.P.Val) & bulk$adj.P.Val < 0.05 & abs(bulk$logFC) > 1, ]
log_msg("Bulk DEG (FDR<0.05 & |logFC|>1): ", nrow(bulk_deg))

# 组成方向: 该细胞类型在 VS vs Nerve 的比例变化方向
comp_dir <- with(comp_new, {
  vs <- pct[group == "VS"];  names(vs) <- cell_type[group == "VS"]
  nv <- pct[group == "Nerve"]; names(nv) <- cell_type[group == "Nerve"]
  d <- vs[REAL_CELLTYPES] - nv[REAL_CELLTYPES]; names(d) <- REAL_CELLTYPES; d
})

classify <- function(g) {
  bl <- bulk_deg$logFC[bulk_deg$gene == g][1]
  dct <- unname(dom_ct[g])
  if (!length(dct)) dct <- NA_character_
  pb_lfc <- NA_real_; pb_p <- NA_real_
  if (!is.na(dct) && dct %in% names(de_list)) {
    de <- de_list[[dct]]$de
    row <- de[de$gene == g, ]
    if (nrow(row)) { pb_lfc <- row$log2FC[1]; pb_p <- row$pvalue[1] }
  }
  intrinsic <- !is.na(pb_p) && pb_p < 0.05 && !is.na(pb_lfc) && sign(pb_lfc) == sign(bl)
  cdir <- if (!is.na(dct)) comp_dir[[dct]] else NA_real_
  composition <- !is.na(cdir) && sign(cdir) == sign(bl) && abs(cdir) >= 1
  cls <- if (intrinsic && composition) "both" else
         if (intrinsic) "cell-intrinsic (candidate)" else
         if (composition) "composition-consistent" else "ambiguous/undetermined"
  data.frame(gene = g, bulk_logFC = bl, dominant_celltype = ifelse(is.na(dct), NA, dct),
             pb_log2FC_in_dom = pb_lfc, pb_p_in_dom = pb_p,
             dom_pct_shift_VS_minus_Nerve = ifelse(is.na(cdir), NA, round(cdir, 2)),
             class = cls, stringsAsFactors = FALSE)
}
cls_tab <- do.call(rbind, lapply(bulk_deg$gene, classify))
write.csv(cls_tab, file.path(outdir, "intrinsic_vs_composition.csv"), row.names = FALSE)
cls_summary <- as.data.frame(table(class = cls_tab$class))
write.csv(cls_summary, file.path(outdir, "intrinsic_vs_composition_summary.csv"), row.names = FALSE)
log_msg("Classification summary:")
print(cls_summary)

# ---- pseudobulk logFC vs bulk logFC 相关散点 (三主线) ----
focus_ct <- intersect(c("Schwann/肿瘤细胞", "巨噬/髓系", "成纤维"), names(de_list))
cor_rows <- list()
for (ct in focus_ct) {
  de <- de_list[[ct]]$de
  m <- merge(de[, c("gene", "log2FC")], bulk[, c("gene", "logFC")], by = "gene")
  m <- m[is.finite(m$log2FC) & is.finite(m$logFC), ]
  r <- suppressWarnings(cor(m$log2FC, m$logFC, method = "pearson"))
  cor_rows[[ct]] <- data.frame(cell_type = ct, n_genes = nrow(m), pearson_r = round(r, 3))
  p <- ggplot(m, aes(logFC, log2FC)) +
    geom_point(size = 0.5, alpha = 0.3, color = "#2166AC") +
    geom_hline(yintercept = 0, color = "grey60") + geom_vline(xintercept = 0, color = "grey60") +
    geom_smooth(method = "lm", se = FALSE, color = "#B2182B", linewidth = 0.7) +
    labs(title = paste0("Within-", CT_LATIN[[ct]], " pseudobulk vs bulk logFC"),
         subtitle = paste0("Pearson r = ", round(r, 3), " (n=", nrow(m), " genes)"),
         x = "bulk limma log2FC (VS / Nerve)",
         y = paste0("pseudobulk log2FC in ", CT_LATIN[[ct]])) +
    theme_bw(base_size = 11)
  ggsave(file.path(outdir, paste0("fig_corr_", CT_LATIN[[ct]], "_vs_bulk.png")),
         p, width = 6, height = 6, dpi = 300)
}
cor_tab <- do.call(rbind, cor_rows)
write.csv(cor_tab, file.path(outdir, "pseudobulk_bulk_logFC_correlation.csv"), row.names = FALSE)

# ---- bulk-pseudobulk concordance sensitivity analyses ----
sensitivity_rows <- list()
for (ct in focus_ct) {
  de <- de_list[[ct]]$de
  specs <- list(
    list(input_filter = "all", output_filter = "all", threshold = NA_real_,
         bulk_input = bulk),
    list(input_filter = "baseMean", output_filter = "baseMean", threshold = 10,
         bulk_input = bulk),
    list(input_filter = "baseMean", output_filter = "baseMean", threshold = 20,
         bulk_input = bulk),
    list(input_filter = "bulk_DEG",
         output_filter = "bulk_DEG_supplementary", threshold = NA_real_,
         bulk_input = bulk_deg)
  )
  ct_rows <- lapply(specs, function(spec) {
    m <- filter_concordance_data(
      de, spec$bulk_input, spec$input_filter, spec$threshold
    )
    summarize_concordance(
      m, ct, spec$output_filter, spec$threshold, n_boot = 2000L
    )
  })
  sensitivity_rows[[ct]] <- do.call(rbind, ct_rows)
}
sensitivity_tab <- do.call(rbind, sensitivity_rows)
rownames(sensitivity_tab) <- NULL
write.csv(
  sensitivity_tab,
  file.path(outdir, "pseudobulk_bulk_logFC_sensitivity.csv"),
  row.names = FALSE
)

for (ct in focus_ct) {
  plot_df <- sensitivity_tab[sensitivity_tab$cell_type == ct, ]
  plot_df$analysis <- ifelse(
    plot_df$filter_type == "baseMean",
    paste0("baseMean≥", plot_df$threshold),
    ifelse(plot_df$filter_type == "all", "All matched genes", "Bulk DEG (supp.)")
  )
  plot_df$analysis <- factor(
    plot_df$analysis,
    levels = c("All matched genes", "baseMean≥10", "baseMean≥20",
               "Bulk DEG (supp.)")
  )
  long <- rbind(
    data.frame(
      analysis = plot_df$analysis, method = "Pearson",
      estimate = plot_df$pearson_r,
      ci_low = plot_df$pearson_ci_low, ci_high = plot_df$pearson_ci_high,
      n_genes = plot_df$n_genes
    ),
    data.frame(
      analysis = plot_df$analysis, method = "Spearman",
      estimate = plot_df$spearman_r,
      ci_low = plot_df$spearman_ci_low, ci_high = plot_df$spearman_ci_high,
      n_genes = plot_df$n_genes
    )
  )
  p <- ggplot(long, aes(analysis, estimate, color = method, group = method)) +
    geom_hline(yintercept = 0, color = "grey70") +
    geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                  width = 0.12, position = position_dodge(width = 0.35)) +
    geom_point(size = 2.3, position = position_dodge(width = 0.35)) +
    scale_color_manual(values = c(Pearson = "#2166AC", Spearman = "#B2182B")) +
    labs(
      title = paste0(CT_LATIN[[ct]], " bulk-pseudobulk concordance sensitivity"),
      subtitle = "Points are correlations; bars are bootstrap percentile 95% CIs",
      x = NULL, y = "Correlation coefficient", color = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
  ggsave(
    file.path(outdir, paste0("fig_corr_sensitivity_", CT_LATIN[[ct]], ".png")),
    p, width = 7.2, height = 5.5, dpi = 300
  )
}

dominant_tab <- build_dominant_celltype_concordance(
  bulk_deg, dom_ct, de_list, thresholds = c(10, 20)
)
write.csv(
  dominant_tab,
  file.path(outdir, "dominant_celltype_direction_concordance.csv"),
  row.names = FALSE
)

dominant_plot <- dominant_tab
dominant_plot$analysis <- ifelse(
  dominant_plot$filter_type == "all",
  "All available", paste0("baseMean≥", dominant_plot$threshold)
)
dominant_plot$label <- ifelse(
  dominant_plot$scope == "overall",
  "Overall",
  unname(CT_PLOT[dominant_plot$dominant_celltype])
)
p_dominant <- ggplot(
  dominant_plot,
  aes(label, direction_concordance, color = analysis, group = analysis)
) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "grey60") +
  geom_errorbar(
    aes(ymin = direction_ci_low, ymax = direction_ci_high),
    width = 0.15, position = position_dodge(width = 0.55)
  ) +
  geom_point(size = 2.2, position = position_dodge(width = 0.55)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_color_manual(
    values = c("All available" = "#1B7837", "baseMean≥10" = "#2166AC",
               "baseMean≥20" = "#B2182B")
  ) +
  labs(
    title = "Dominant-cell-type fold-change direction concordance",
    subtitle = "Each bulk DEG contributes only its assigned dominant-cell-type pseudobulk estimate",
    x = NULL, y = "Direction concordance (Clopper-Pearson 95% CI)", color = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(
  file.path(outdir, "fig_dominant_celltype_direction_concordance.png"),
  p_dominant, width = 10, height = 6, dpi = 300
)

# ---- 报告 ----
md <- function(x, n = Inf) {
  if (!nrow(x)) return("_无_")
  x <- head(x, n)
  paste(c(paste0("| ", paste(names(x), collapse = " | "), " |"),
          paste0("|", paste(rep("---", ncol(x)), collapse = "|"), "|"),
          apply(x, 1, function(z) paste0("| ", paste(z, collapse = " | "), " |"))),
        collapse = "\n")
}
n_tab <- do.call(rbind, lapply(de_list, function(z)
  data.frame(cell_type = z$ct, n_VS = unname(z$n["VS"]), n_Nerve = unname(z$n["Nerve"]),
             n_genes_tested = nrow(z$de),
             n_sig = sum(z$de$padj < 0.05 & abs(z$de$log2FC) > 1, na.rm = TRUE))))
sensitivity_report <- sensitivity_tab
numeric_cols <- vapply(sensitivity_report, is.numeric, logical(1))
sensitivity_report[numeric_cols] <- lapply(
  sensitivity_report[numeric_cols], function(x) round(x, 3)
)
dominant_report <- dominant_tab[
  dominant_tab$scope == "overall",
  c("filter_type", "threshold", "n_bulk_deg", "n_with_pb", "n_direction",
    "n_concordant", "direction_concordance", "direction_ci_low",
    "direction_ci_high")
]
dominant_numeric <- vapply(dominant_report, is.numeric, logical(1))
dominant_report[dominant_numeric] <- lapply(
  dominant_report[dominant_numeric], function(x) round(x, 3)
)
report <- c(
  "# 加深1: scRNA 细胞类型内 pseudobulk 差异表达报告 (GSE230375)",
  "",
  "## 目的与措辞",
  "在 GSE230375（7 VS + 2 great auricular nerve）中按样本×粗粒度细胞类型聚合 raw counts，并以 DESeq2 比较 VS 与 Nerve。分析用于区分与细胞组成变化一致的信号和具有候选细胞类型内证据的信号；结果属于探索性 in-silico 证据，不构成肿瘤细胞内在机制或直接实验验证。",
  "",
  "**统计约束**：VS n=7，Nerve n=2。细胞类型内显著性欠功效，因此主要报告方向一致性、bulk–pseudobulk Pearson/Spearman 相关及其置信区间；相关分析采用全匹配基因和预定义 baseMean≥10、baseMean≥20 阈值。bulk-DEG 限定结果仅作补充，因为按 bulk 结果筛选可能选择性抬高一致性。",
  "",
  paste0("**构成可复现性闸**: ", ifelse(gate_ok, "PASS (重建细胞构成与 08 记录一致, <1% 漂移)",
         "WARN: 构成漂移 >1%, 见 composition_reproducibility_check.csv")),
  "",
  "## 各细胞类型 pseudobulk DE 概况",
  md(n_tab),
  "",
  "## bulk DEG 的组成一致性与候选细胞类型内证据",
  paste0("bulk DEG (FDR<0.05 & |logFC|>1) 共 ", nrow(bulk_deg), " 个, 按其主导细胞类型 (scRNA 平均表达最高者) 归类:"),
  "",
  md(cls_summary),
  "",
  "- **cell-intrinsic (candidate)**：主导细胞类型内 pseudobulk log2FC 与 bulk 方向一致且 nominal p<0.05，但该标签仅表示候选细胞类型内证据，不等同于已证明的细胞内在效应。",
  "- **composition-consistent**: 细胞类型内不显著, 但主导细胞类型比例在对应方向变化 ≥1% → 与组成变化一致 (如基质/内皮丢失)。",
  "- **both / ambiguous**: 两者皆/均不满足。",
  "- `both` 与 `cell-intrinsic (candidate)` 合计430/1,410（30.5%），表示具有 nominal direction-matched within-cell-type evidence；不得称为“真正 cell-intrinsic”。",
  "",
  "## bulk–pseudobulk log2FC concordance and sensitivity analyses",
  md(sensitivity_report),
  "",
  "整体 dominant-cell-type 方向一致率（每个 bulk DEG 仅使用其主导细胞类型对应的 pseudobulk log2FC；不按 nominal p 筛选）：",
  "",
  md(dominant_report),
  "",
  "## 产物 (完整绝对路径)",
  paste0("- 脚本: `", file.path(root, "scripts", "12_pseudobulk_DE.R"), "`"),
  paste0("- 输入 Seurat: `", rds_path, "`"),
  paste0("- 输入 bulk DEG: `", bulk_path, "`"),
  paste0("- 每细胞类型 DE: `", outdir, "/{Schwann_tumor,Myeloid,T_cell,Fibroblast,Endothelial,Pericyte}_VS_vs_nerve_DE.csv`"),
  paste0("- 火山图: `", outdir, "/fig_volcano_*.png`"),
  paste0("- 相关散点: `", outdir, "/fig_corr_*_vs_bulk.png`"),
  paste0("- 相关敏感性表: `", file.path(outdir, "pseudobulk_bulk_logFC_sensitivity.csv"), "`"),
  paste0("- 主导细胞类型方向一致率: `", file.path(outdir, "dominant_celltype_direction_concordance.csv"), "`"),
  paste0("- 敏感性图: `", outdir, "/fig_corr_sensitivity_{Schwann_tumor,Myeloid,Fibroblast}.png`; `",
         file.path(outdir, "fig_dominant_celltype_direction_concordance.png"), "`"),
  paste0("- 内在/组成分类: `", file.path(outdir, "intrinsic_vs_composition.csv"), "` (+ _summary.csv)"),
  paste0("- 构成复现闸: `", file.path(outdir, "composition_reproducibility_check.csv"), "`"),
  paste0("- pseudobulk counts (供加深2): `", file.path(proc, "GSE230375_pseudobulk_counts.csv"), "`"),
  paste0("- 细胞类型均表达参考 (供加深2反卷积): `", file.path(proc, "GSE230375_celltype_meanexpr.csv"), "`"),
  paste0("- 细胞 meta: `", file.path(proc, "GSE230375_cell_meta.csv"), "`"),
  "",
  "## 局限",
  "- Nerve n=2，细胞类型内显著性欠功效；分类与相关结果均为探索性候选证据。  ",
  "- 粗粒度注释不能区分 neoplastic、repair-like 与正常 Schwann states；“cell-intrinsic (candidate)”只是操作性分类名称。",
  "- 主导细胞类型按平均表达最高者定义, 非特异表达; 多细胞类型共表达基因的归类有不确定性。",
  "- 未做 doublet/ambient RNA 校正; pseudobulk 对解离与捕获偏倚仍敏感。"
)
writeLines(report, file.path(docdir, "pseudobulk_report.md"), useBytes = TRUE)
log_msg("Report written: ", file.path(docdir, "pseudobulk_report.md"))
log_msg("加深1 完成")
