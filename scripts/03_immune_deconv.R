#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(GSVA)
  library(ggplot2)
  library(pheatmap)
})

options(stringsAsFactors = FALSE)
set.seed(20260621)

ROOT <- "."
OUTDIR <- file.path(ROOT, "results/immune")
REPORT_PATH <- file.path(ROOT, "docs/03_immune_report.md")

get_marker_sets <- function() {
  # Bindea/Charoentong-style immune marker sets. The 28 immune populations are
  # supplemented with fibroblast and endothelial stromal marker sets.
  list(
    "B cells" = c("CD19", "MS4A1", "CD79A", "CD79B", "CD37", "CD74",
                  "HLA-DRA", "CD22", "BANK1", "BLK"),
    "Memory B cells" = c("CD27", "TNFRSF13B", "GPR183", "CD44", "CD82",
                         "AIM2", "BANK1", "CD79A"),
    "Plasma cells" = c("MZB1", "SDC1", "CD38", "JCHAIN", "DERL3", "XBP1",
                       "TNFRSF17", "IGHG1", "IGHG3"),
    "T cells" = c("CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2",
                  "LCK", "CD247"),
    "CD8 T cells" = c("CD8A", "CD8B", "GZMK", "CCL5", "NKG7", "LCK",
                      "TRAC", "CD3D"),
    "CD4 T cells" = c("CD4", "IL7R", "LTB", "MAL", "NOSIP", "TRAT1",
                      "CD3D", "CD3E"),
    "T helper cells" = c("CD4", "IL7R", "MAL", "LTB", "ICOS", "CD40LG",
                         "TRAT1", "CD3D"),
    "Th1 cells" = c("TBX21", "STAT4", "IL12RB2", "IFNG", "CXCR3", "IL18R1",
                    "CCL5", "TNF"),
    "Th2 cells" = c("GATA3", "STAT6", "IL4R", "CCR4", "IL17RB", "GPR44",
                    "ICOS", "IL13RA1"),
    "Th17 cells" = c("RORC", "CCR6", "KLRB1", "IL23R", "RORA", "AHR",
                     "IL7R", "CCL20"),
    "Treg cells" = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT", "CCR8",
                     "TNFRSF18", "ENTPD1"),
    "T follicular helper cells" = c("CXCR5", "PDCD1", "ICOS", "BCL6",
                                    "IL21R", "SH2D1A", "TOX", "MAF"),
    "Gamma delta T cells" = c("TRDC", "TRGC1", "TRGC2", "CD3D", "CD3E",
                              "KLRD1", "NKG7", "CCL5"),
    "Cytotoxic cells" = c("NKG7", "GNLY", "GZMB", "GZMH", "PRF1", "CTSW",
                          "FGFBP2", "KLRD1"),
    "NK cells" = c("KLRD1", "FCGR3A", "NKG7", "GNLY", "PRF1", "TRAC",
                   "TYROBP", "XCL1"),
    "CD56bright NK cells" = c("NCAM1", "XCL1", "XCL2", "KLRC1", "KLRD1",
                              "TRAC", "IL7R", "GZMK"),
    "CD56dim NK cells" = c("FCGR3A", "FGFBP2", "CX3CR1", "PRF1", "GZMB",
                           "KLRD1", "NKG7", "GNLY"),
    "Dendritic cells" = c("FCER1A", "CD1C", "CST3", "CLEC10A", "HLA-DPA1",
                          "HLA-DPB1", "HLA-DRA", "ITGAX"),
    "Activated dendritic cells" = c("CD83", "CD86", "CCR7", "LAMP3",
                                    "HLA-DRA", "HLA-DPA1", "IL15", "MARCKSL1"),
    "Immature dendritic cells" = c("CD1A", "CD1C", "FCER1A", "CLEC10A",
                                   "CST3", "HLA-DRA", "IRF4", "ITGAX"),
    "Plasmacytoid dendritic cells" = c("GZMB", "JCHAIN", "TCF4", "CLEC4C",
                                      "IL3RA", "IRF7", "TSPAN13", "PACSIN1",
                                      "LILRB4", "IGJ"),
    "Macrophages" = c("CD68", "CD163", "CSF1R", "C1QA", "C1QB", "C1QC",
                      "MSR1", "CTSD", "FCER1G", "TYROBP"),
    "Monocytes" = c("S100A8", "S100A9", "CTSS", "FCN1", "VCAN", "LYZ",
                    "CTSD", "FCGR3A"),
    "MDSC" = c("S100A8", "S100A9", "ARG1", "IL4R", "CYBB", "ITGAM",
               "FCGR3B", "CEACAM8"),
    "Neutrophils" = c("FCGR3B", "CSF3R", "S100A8", "S100A9", "CEACAM8",
                      "FPR1", "CXCR2", "MMP8"),
    "Eosinophils" = c("CLC", "SIGLEC8", "IL5RA", "CCR3", "PRG2", "RNASE2",
                      "RNASE3", "ALOX15"),
    "Mast cells" = c("TPSAB1", "TPSB2", "KIT", "MS4A2", "CPA3", "HDC",
                     "HPGDS", "GATA2"),
    "Basophils" = c("HDC", "MS4A2", "FCER1A", "CCR3", "IL3RA", "CLC",
                    "GATA2", "ENPP3"),
    "Fibroblasts" = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "COL6A1",
                      "PDGFRA", "FAP", "THY1", "COL5A1"),
    "Endothelial cells" = c("PECAM1", "VWF", "KDR", "EMCN", "ENG", "RAMP2",
                            "PLVAP", "ESAM", "CDH5", "EGFL7")
  )
}

validate_and_align_inputs <- function(expr, meta) {
  stopifnot(
    is.matrix(expr),
    !is.null(rownames(expr)),
    !is.null(colnames(expr)),
    all(c("sample", "group", "batch") %in% names(meta)),
    !anyDuplicated(rownames(expr)),
    !anyDuplicated(colnames(expr)),
    !anyDuplicated(meta$sample)
  )
  meta <- as.data.table(meta)
  meta <- meta[match(colnames(expr), sample)]
  stopifnot(
    !anyNA(meta$sample),
    identical(meta$sample, colnames(expr)),
    all(meta$group %in% c("Control", "Tumor")),
    !anyNA(meta$batch)
  )
  meta[, group := factor(group, levels = c("Control", "Tumor"))]
  meta[, batch := factor(batch)]
  list(expr = expr, meta = meta)
}

compare_scores <- function(scores, meta) {
  meta <- as.data.table(meta)
  meta <- meta[match(colnames(scores), sample)]
  rows <- lapply(seq_len(nrow(scores)), function(i) {
    control <- as.numeric(scores[i, meta$group == "Control"])
    vs <- as.numeric(scores[i, meta$group == "Tumor"])
    wt <- wilcox.test(vs, control, exact = FALSE)
    data.table(
      cell_type = rownames(scores)[i],
      n_control = length(control),
      n_VS = length(vs),
      median_control = median(control, na.rm = TRUE),
      median_VS = median(vs, na.rm = TRUE),
      delta_median = median(vs, na.rm = TRUE) - median(control, na.rm = TRUE),
      p_value = wt$p.value
    )
  })
  ans <- rbindlist(rows)
  ans[, p_adj := p.adjust(p_value, method = "BH")]
  ans[, direction := fifelse(
    delta_median > 0, "Higher in VS",
    fifelse(delta_median < 0, "Lower in VS", "No median difference")
  )]
  ans[, abs_delta_median := abs(delta_median)]
  setorder(ans, p_adj, -abs_delta_median)
  ans[, abs_delta_median := NULL]
  ans
}

format_num <- function(x, digits = 3L) {
  ifelse(
    is.na(x), "NA",
    ifelse(abs(x) < 0.001, format(x, scientific = TRUE, digits = digits),
           format(round(x, digits), nsmall = digits, trim = TRUE))
  )
}

run_analysis <- function() {
  dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(REPORT_PATH), recursive = TRUE, showWarnings = FALSE)

  expr_dt <- fread(file.path(ROOT, "data/processed/GSE39645_expr_gene.csv"))
  meta <- fread(file.path(ROOT, "data/processed/GSE39645_meta.csv"))
  signature <- fread(file.path(ROOT, "results/signature/signature_genes.csv"))

  stopifnot("gene" %in% names(expr_dt), "gene" %in% names(signature))
  expr <- as.matrix(expr_dt[, -1])
  storage.mode(expr) <- "double"
  rownames(expr) <- expr_dt$gene
  validated <- validate_and_align_inputs(expr, meta)
  expr <- validated$expr
  meta <- validated$meta

  group_design <- model.matrix(~ group, data = meta)
  expr_adj <- removeBatchEffect(
    expr,
    batch = meta$batch,
    design = group_design
  )

  marker_sets <- lapply(get_marker_sets(), unique)
  marker_coverage <- rbindlist(lapply(names(marker_sets), function(cell_type) {
    data.table(
      cell_type = cell_type,
      marker_count_defined = length(marker_sets[[cell_type]]),
      marker_count_detected = sum(marker_sets[[cell_type]] %in% rownames(expr_adj)),
      detected_markers = paste(
        intersect(marker_sets[[cell_type]], rownames(expr_adj)),
        collapse = ";"
      )
    )
  }))
  fwrite(marker_coverage, file.path(OUTDIR, "marker_coverage.csv"))

  usable_sets <- lapply(marker_sets, intersect, y = rownames(expr_adj))
  usable_sets <- usable_sets[lengths(usable_sets) >= 5L]
  stopifnot(length(usable_sets) >= 28L)

  # GSVA >= 2.0 parameter API; this is ssGSEA with the historical
  # gsva(..., method = "ssgsea") behavior.
  ssgsea_param <- ssgseaParam(
    exprData = expr_adj,
    geneSets = usable_sets,
    minSize = 5L,
    maxSize = Inf,
    normalize = TRUE,
    checkNA = "auto"
  )
  ssgsea_scores <- gsva(ssgsea_param, verbose = FALSE)
  ssgsea_scores <- as.matrix(ssgsea_scores)
  score_out <- data.table(cell_type = rownames(ssgsea_scores), ssgsea_scores)
  fwrite(score_out, file.path(OUTDIR, "ssgsea_scores.csv"))

  diff <- compare_scores(ssgsea_scores, meta)
  fwrite(diff, file.path(OUTDIR, "immune_diff_VS_vs_control.csv"))
  significant_cells <- diff[p_adj < 0.05, cell_type]

  plot_cells <- significant_cells
  plot_subtitle <- "BH-adjusted p < 0.05"
  if (!length(plot_cells)) {
    plot_cells <- head(diff$cell_type, 6L)
    plot_subtitle <- "No BH-significant populations; six lowest adjusted p-values shown"
  }
  score_long <- melt(
    score_out[cell_type %in% plot_cells],
    id.vars = "cell_type",
    variable.name = "sample",
    value.name = "ssGSEA_score"
  )
  score_long <- merge(
    score_long,
    meta[, .(sample, group = fifelse(group == "Tumor", "VS", "Control"))],
    by = "sample"
  )
  score_long[, cell_type := factor(cell_type, levels = plot_cells)]
  p_box <- ggplot(score_long, aes(x = group, y = ssGSEA_score, fill = group)) +
    geom_boxplot(width = 0.65, outlier.shape = NA, alpha = 0.82) +
    geom_jitter(width = 0.14, size = 1.25, alpha = 0.72) +
    facet_wrap(~ cell_type, scales = "free_y", ncol = 4) +
    scale_fill_manual(values = c(Control = "#4C78A8", VS = "#E45756")) +
    labs(
      title = "GSE39645 batch-adjusted ssGSEA scores",
      subtitle = plot_subtitle,
      x = NULL,
      y = "Normalized ssGSEA enrichment score"
    ) +
    theme_bw(base_size = 11) +
    theme(
      legend.position = "none",
      strip.background = element_rect(fill = "grey94"),
      axis.text.x = element_text(angle = 20, hjust = 1)
    )
  ggsave(
    file.path(OUTDIR, "fig_immune_boxplots.png"),
    p_box,
    width = 13,
    height = max(5.5, 3.1 * ceiling(length(plot_cells) / 4)),
    dpi = 180
  )

  log_lines <- c(
    sprintf("Analysis date: %s", Sys.Date()),
    sprintf("GSVA version: %s", as.character(packageVersion("GSVA"))),
    sprintf("Usable ssGSEA gene sets: %d", length(usable_sets))
  )
  mcp_path <- file.path(OUTDIR, "mcpcounter_scores.csv")
  if (requireNamespace("MCPcounter", quietly = TRUE)) {
    mcp_scores <- MCPcounter::MCPcounter.estimate(
      expr_adj,
      featuresType = "HUGO_symbols",
      probesets = rownames(expr_adj)
    )
    fwrite(
      data.table(cell_type = rownames(mcp_scores), mcp_scores),
      mcp_path
    )
    log_lines <- c(log_lines, sprintf("MCPcounter: completed (%s)", mcp_path))
  } else {
    log_lines <- c(
      log_lines,
      "MCPcounter: skipped because the MCPcounter package was not installed."
    )
  }

  sig_genes <- intersect(unique(signature$gene), rownames(expr_adj))
  if (!length(sig_genes)) {
    stop("No signature genes were present in the discovery expression matrix.")
  }
  cor_cells <- significant_cells
  if (length(cor_cells)) {
    cor_rows <- rbindlist(lapply(sig_genes, function(gene) {
      rbindlist(lapply(cor_cells, function(cell_type) {
        ct <- cor.test(
          as.numeric(expr_adj[gene, ]),
          as.numeric(ssgsea_scores[cell_type, ]),
          method = "spearman",
          exact = FALSE
        )
        data.table(
          gene = gene,
          cell_type = cell_type,
          spearman_rho = unname(ct$estimate),
          p_value = ct$p.value
        )
      }))
    }))
    cor_rows[, p_adj := p.adjust(p_value, method = "BH")]
    setorder(cor_rows, gene, p_adj)
  } else {
    cor_rows <- data.table(
      gene = character(),
      cell_type = character(),
      spearman_rho = numeric(),
      p_value = numeric(),
      p_adj = numeric()
    )
  }
  fwrite(cor_rows, file.path(OUTDIR, "signature_immune_cor.csv"))

  heatmap_path <- file.path(OUTDIR, "fig_signature_immune_cor.png")
  if (nrow(cor_rows)) {
    cor_mat <- dcast(
      cor_rows,
      gene ~ cell_type,
      value.var = "spearman_rho"
    )
    mat <- as.matrix(cor_mat[, -1])
    rownames(mat) <- cor_mat$gene
    p_mat <- dcast(cor_rows, gene ~ cell_type, value.var = "p_adj")
    p_values <- as.matrix(p_mat[, -1])
    rownames(p_values) <- p_mat$gene
    stars <- matrix(
      ifelse(
        p_values < 0.001, "***",
        ifelse(p_values < 0.01, "**", ifelse(p_values < 0.05, "*", ""))
      ),
      nrow = nrow(p_values),
      dimnames = dimnames(p_values)
    )
    png(
      heatmap_path,
      width = max(1200, 210 * ncol(mat)),
      height = max(900, 90 * nrow(mat)),
      res = 160
    )
    pheatmap(
      mat,
      color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
      breaks = seq(-1, 1, length.out = 102),
      cluster_rows = nrow(mat) > 1,
      cluster_cols = ncol(mat) > 1,
      display_numbers = stars,
      number_color = "black",
      fontsize_number = 10,
      main = "Signature gene expression vs significant ssGSEA scores\nSpearman correlation; BH-adjusted significance"
    )
    dev.off()
  } else {
    png(heatmap_path, width = 1200, height = 700, res = 160)
    plot.new()
    text(
      0.5, 0.55,
      "No BH-significant ssGSEA population was available\nfor signature–infiltration correlation analysis.",
      cex = 1.2
    )
    text(0.5, 0.42, "The correlation CSV contains headers only.", cex = 0.95)
    dev.off()
  }

  writeLines(log_lines, file.path(OUTDIR, "analysis_log.txt"))

  higher <- diff[p_adj < 0.05 & delta_median > 0]
  lower <- diff[p_adj < 0.05 & delta_median < 0]
  report_lines <- c(
    "# GSE39645 免疫与基质浸润评分报告",
    "",
    "## 方法",
    "",
    paste0(
      "本分析读取 GSE39645 基因级 log2 表达矩阵，并以组别设计矩阵保留 ",
      "VS（元数据中标记为 Tumor）与 Control 的组间差异，同时使用 ",
      "`limma::removeBatchEffect` 校正 batch。随后基于脚本内硬编码的 ",
      "28 类免疫细胞标志基因集及 2 类基质标志基因集进行 ssGSEA。组间比较采用",
      "双侧 Wilcoxon 秩和检验，并使用 Benjamini–Hochberg 方法校正多重检验。"
    ),
    "",
    "## 结果",
    ""
  )
  if (nrow(higher)) {
    higher_text <- paste(
      sprintf(
        "%s（VS 中位数 %s，Control 中位数 %s，Δ中位数 %s，BH 校正 p=%s）",
        higher$cell_type,
        format_num(higher$median_VS),
        format_num(higher$median_control),
        format_num(higher$delta_median),
        format_num(higher$p_adj)
      ),
      collapse = "；"
    )
    report_lines <- c(
      report_lines,
      paste0("在 BH 校正 p<0.05 的标准下，VS 组评分较高的群体包括：", higher_text, "。"),
      ""
    )
  } else {
    report_lines <- c(
      report_lines,
      "在 BH 校正 p<0.05 的标准下，未观察到 VS 组评分显著升高的群体。",
      ""
    )
  }
  if (nrow(lower)) {
    lower_text <- paste(
      sprintf(
        "%s（VS 中位数 %s，Control 中位数 %s，Δ中位数 %s，BH 校正 p=%s）",
        lower$cell_type,
        format_num(lower$median_VS),
        format_num(lower$median_control),
        format_num(lower$delta_median),
        format_num(lower$p_adj)
      ),
      collapse = "；"
    )
    report_lines <- c(
      report_lines,
      paste0("VS 组评分较低的群体包括：", lower_text, "。"),
      ""
    )
  } else {
    report_lines <- c(
      report_lines,
      "在 BH 校正 p<0.05 的标准下，未观察到 VS 组评分显著降低的群体。",
      ""
    )
  }

  if (nrow(cor_rows)) {
    cor_sig <- cor_rows[p_adj < 0.05]
    if (nrow(cor_sig)) {
      top_cor <- cor_sig[order(p_adj, -abs(spearman_rho))]
      top_cor <- head(top_cor, 20L)
      cor_text <- paste(
        sprintf(
          "%s–%s：ρ=%s，BH 校正 p=%s（%s相关）",
          top_cor$gene,
          top_cor$cell_type,
          format_num(top_cor$spearman_rho),
          format_num(top_cor$p_adj),
          ifelse(top_cor$spearman_rho > 0, "正", "负")
        ),
        collapse = "；"
      )
      report_lines <- c(
        report_lines,
        paste0(
          "在 signature 基因与上述差异群体评分的 Spearman 分析中，",
          "BH 校正 p<0.05 的相关包括：", cor_text, "。"
        ),
        ""
      )
    } else {
      report_lines <- c(
        report_lines,
        paste0(
          "signature 基因与差异群体评分之间未观察到 BH 校正 p<0.05 的",
          "Spearman 相关。"
        ),
        ""
      )
    }
  } else {
    report_lines <- c(
      report_lines,
      paste0(
        "由于不存在 BH 校正 p<0.05 的差异群体，未进一步计算 signature ",
        "基因与群体评分的相关。"
      ),
      ""
    )
  }

  report_lines <- c(
    report_lines,
    "## 解释边界",
    "",
    paste0(
      "上述结果来自 bulk 转录组的基因集富集与计算反卷积推断，反映标志基因",
      "表达模式的相对变化，不等同于组织中细胞数量或比例的直接测量。batch 与",
      "组别存在一定不平衡，且 Control 包含不同正常神经及 Schwann 细胞来源，",
      "因此这些结果仅作为次级支撑证据。signature–评分相关是在全队列中计算，",
      "可能部分反映两者共同的组间差异，不能据此推断细胞内调控或直接作用。",
      "MCPcounter 是否执行见分析日志。"
    ),
    "",
    "## 产物路径",
    "",
    paste0("ssGSEA 分数矩阵：`", file.path(OUTDIR, "ssgsea_scores.csv"), "`"),
    "",
    paste0("标志基因覆盖：`", file.path(OUTDIR, "marker_coverage.csv"), "`"),
    "",
    paste0("VS vs Control 差异结果：`", file.path(OUTDIR, "immune_diff_VS_vs_control.csv"), "`"),
    "",
    paste0("差异群体箱线图：`", file.path(OUTDIR, "fig_immune_boxplots.png"), "`"),
    "",
    paste0("signature–群体相关：`", file.path(OUTDIR, "signature_immune_cor.csv"), "`"),
    "",
    paste0("signature–群体相关热图：`", heatmap_path, "`"),
    "",
    paste0("运行日志：`", file.path(OUTDIR, "analysis_log.txt"), "`"),
    "",
    paste0("分析脚本：`", file.path(ROOT, "scripts/03_immune_deconv.R"), "`"),
    "",
    paste0("契约测试：`", file.path(ROOT, "scripts/test_03_immune_deconv.R"), "`")
  )
  if (file.exists(mcp_path)) {
    report_lines <- c(
      report_lines,
      "",
      paste0("MCPcounter 分数：`", mcp_path, "`")
    )
  }
  writeLines(report_lines, REPORT_PATH)

  cat(sprintf("ssGSEA populations: %d\n", nrow(ssgsea_scores)))
  cat(sprintf("BH-significant populations: %d\n", length(significant_cells)))
  cat(sprintf("Signature genes detected: %d/%d\n", length(sig_genes), nrow(signature)))
  cat(sprintf("Report: %s\n", REPORT_PATH))
}

if (sys.nframe() == 0L) {
  run_analysis()
}
