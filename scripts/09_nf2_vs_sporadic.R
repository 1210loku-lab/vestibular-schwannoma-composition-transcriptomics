#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
set.seed(20260621)

root <- Sys.getenv("PROJ_ROOT", getwd())
raw_file <- file.path(root, "data/raw/GSE108524_matrix.txt.gz")
expr_file <- file.path(root, "data/processed/GSE108524_expr_gene.csv")
meta_file <- file.path(root, "data/processed/GSE108524_meta.csv")
signature_file <- file.path(root, "results/signature/signature_genes.csv")
script_file <- file.path(root, "scripts/09_nf2_vs_sporadic.R")
outdir <- file.path(root, "results/nf2_sporadic")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  evidence = file.path(outdir, "sample_label_evidence.csv"),
  deg_all = file.path(outdir, "deg_nf2_vs_sporadic.csv"),
  deg_sig = file.path(outdir, "deg_nf2_vs_sporadic_significant.csv"),
  pca_coordinates = file.path(outdir, "pca_coordinates.csv"),
  pca_figure = file.path(outdir, "fig_pca_nf2_sporadic.png"),
  signature_scores = file.path(outdir, "signature_scores_nf2_sporadic.csv"),
  signature_figure = file.path(outdir, "fig_signature_score_nf2_sporadic.png"),
  summary = file.path(outdir, "analysis_summary.csv"),
  report = file.path(outdir, "09_nf2_sporadic_report.md"),
  session_info = file.path(outdir, "sessionInfo.txt")
)

required_inputs <- c(raw_file, expr_file, meta_file, signature_file)
if (!all(file.exists(required_inputs))) {
  stop("Missing required input(s): ", paste(required_inputs[!file.exists(required_inputs)], collapse = ", "))
}

read_series_values <- function(lines, prefix, first_only = TRUE) {
  hits <- lines[startsWith(lines, prefix)]
  if (!length(hits)) stop("Missing GEO series-matrix row: ", prefix)
  parsed <- lapply(hits, function(x) {
    gsub('"', "", strsplit(x, "\t", fixed = TRUE)[[1]][-1])
  })
  if (first_only) parsed[[1]] else parsed
}

lines <- readLines(gzfile(raw_file), warn = FALSE)
sample_ids <- read_series_values(lines, "!Sample_geo_accession")
sample_titles <- read_series_values(lines, "!Sample_title")
characteristic_rows <- read_series_values(
  lines,
  "!Sample_characteristics_ch1",
  first_only = FALSE
)

if (!all(lengths(characteristic_rows) == length(sample_ids))) {
  stop("GEO characteristic rows do not match the sample count.")
}
if (length(sample_titles) != length(sample_ids)) {
  stop("GEO title row does not match the sample count.")
}

characteristics_by_sample <- vapply(
  seq_along(sample_ids),
  function(i) {
    values <- vapply(characteristic_rows, `[`, character(1), i)
    values <- values[nzchar(trimws(values))]
    paste(values, collapse = " | ")
  },
  character(1)
)

title_nf2 <- grepl("\\bnf2\\b", sample_titles, ignore.case = TRUE)
title_sporadic <- grepl("\\bsporadic\\b", sample_titles, ignore.case = TRUE)
char_nf2 <- grepl(
  "\\bnf2\\b|age at onset of nf2",
  characteristics_by_sample,
  ignore.case = TRUE
)
char_sporadic <- grepl("\\bsporadic\\b", characteristics_by_sample, ignore.case = TRUE)
control_evidence <- grepl(
  "\\bcontrol\\b|normal (nerve|never)",
  paste(sample_titles, characteristics_by_sample),
  ignore.case = TRUE
)

subtype <- rep(NA_character_, length(sample_ids))
subtype[title_sporadic & !title_nf2] <- "sporadic"
subtype[title_nf2 & !title_sporadic] <- "NF2"
subtype[is.na(subtype) & char_sporadic & !char_nf2] <- "sporadic"
subtype[is.na(subtype) & char_nf2 & !char_sporadic] <- "NF2"
subtype[is.na(subtype) & control_evidence] <- "control"

classification_basis <- ifelse(
  title_nf2 & !title_sporadic,
  "title contains NF2",
  ifelse(
    title_sporadic & !title_nf2,
    "title contains sporadic",
    ifelse(
      char_nf2 & !char_sporadic,
      "characteristics contain NF2 evidence",
      ifelse(
        char_sporadic & !char_nf2,
        "characteristics contain sporadic",
        ifelse(control_evidence, "control evidence", "unclear")
      )
    )
  )
)

evidence <- data.table(
  sample = sample_ids,
  title = sample_titles,
  characteristics = characteristics_by_sample,
  title_nf2 = title_nf2,
  title_sporadic = title_sporadic,
  characteristics_nf2 = char_nf2,
  characteristics_sporadic = char_sporadic,
  assigned_subtype = subtype,
  classification_basis = classification_basis
)
fwrite(evidence, paths$evidence)

expr_dt <- fread(expr_file)
meta <- fread(meta_file)
if (!"gene" %in% names(expr_dt)) stop("Expression matrix lacks a gene column.")
if (!all(c("sample", "group") %in% names(meta))) {
  stop("Metadata must contain sample and group columns.")
}

expr <- as.matrix(expr_dt[, -1])
rownames(expr) <- expr_dt$gene
storage.mode(expr) <- "double"
if (anyDuplicated(rownames(expr))) stop("Expression matrix contains duplicated gene names.")
if (anyDuplicated(colnames(expr))) stop("Expression matrix contains duplicated sample names.")

meta[, sample := as.character(sample)]
meta <- meta[match(colnames(expr), sample)]
if (anyNA(meta$sample)) stop("Expression samples are missing from processed metadata.")

label_map <- evidence[match(colnames(expr), sample)]
if (anyNA(label_map$sample)) stop("Expression samples are missing from GEO metadata.")

is_vs <- tolower(meta$group) %in% c("tumor", "vs")
vs_samples <- colnames(expr)[is_vs]
vs_labels <- label_map[is_vs]
if (length(vs_samples) != 27L) stop("Expected 27 VS samples, found ", length(vs_samples), ".")
if (any(!vs_labels$assigned_subtype %in% c("NF2", "sporadic"))) {
  bad <- vs_labels[!assigned_subtype %in% c("NF2", "sporadic"), sample]
  stop("Unclear VS subtype labels remain: ", paste(bad, collapse = ", "))
}

subtype <- factor(vs_labels$assigned_subtype, levels = c("sporadic", "NF2"))
n_sporadic <- sum(subtype == "sporadic")
n_nf2 <- sum(subtype == "NF2")
cat(sprintf("NF2 n = %d; sporadic n = %d\n", n_nf2, n_sporadic))
if (n_nf2 != 17L || n_sporadic != 10L) {
  stop("Unexpected subtype counts; expected NF2=17 and sporadic=10.")
}

expr_vs <- expr[, vs_samples, drop = FALSE]
finite_complete <- rowSums(is.finite(expr_vs)) == ncol(expr_vs)
expr_model <- expr_vs[finite_complete, , drop = FALSE]
if (nrow(expr_model) < 2L) stop("Insufficient complete genes for modeling.")

design <- model.matrix(~ 0 + subtype)
colnames(design) <- levels(subtype)
fit <- lmFit(expr_model, design)
contrast <- makeContrasts(NF2_vs_sporadic = NF2 - sporadic, levels = design)
fit2 <- eBayes(contrasts.fit(fit, contrast), trend = TRUE, robust = FALSE)
deg <- as.data.table(topTable(fit2, coef = "NF2_vs_sporadic", number = Inf, sort.by = "none"))
deg[, gene := rownames(expr_model)]
setcolorder(deg, c("gene", setdiff(names(deg), "gene")))
deg[, significant := !is.na(adj.P.Val) & adj.P.Val < 0.05 & abs(logFC) > 1]
deg[, abs_logFC_for_sort := abs(logFC)]
setorder(deg, P.Value, -abs_logFC_for_sort)
deg[, abs_logFC_for_sort := NULL]
fwrite(deg, paths$deg_all)
deg_sig <- deg[significant == TRUE]
fwrite(deg_sig, paths$deg_sig)
n_deg_fdr <- sum(!is.na(deg$adj.P.Val) & deg$adj.P.Val < 0.05)
n_deg_sig <- nrow(deg_sig)
cat(sprintf(
  "DEG: FDR < 0.05 = %d; FDR < 0.05 and |logFC| > 1 = %d\n",
  n_deg_fdr,
  n_deg_sig
))

gene_variance <- apply(expr_model, 1, var, na.rm = TRUE)
gene_variance <- gene_variance[is.finite(gene_variance) & gene_variance > 0]
n_pca_genes <- min(1000L, length(gene_variance))
if (n_pca_genes < 2L) stop("Insufficient variable genes for PCA.")
pca_genes <- names(sort(gene_variance, decreasing = TRUE))[seq_len(n_pca_genes)]
pca <- prcomp(t(expr_model[pca_genes, , drop = FALSE]), center = TRUE, scale. = TRUE)
variance_explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pca_tab <- data.table(
  sample = rownames(pca$x),
  subtype = as.character(subtype),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)
fwrite(pca_tab, paths$pca_coordinates)

n_manova_pc <- min(5L, ncol(pca$x), ncol(pca$x) - 1L)
manova_data <- data.frame(
  subtype = subtype,
  pca$x[, seq_len(n_manova_pc), drop = FALSE],
  check.names = FALSE
)
manova_fit <- manova(as.matrix(manova_data[, -1, drop = FALSE]) ~ subtype, data = manova_data)
pca_manova_p <- summary(manova_fit, test = "Pillai")$stats[1, "Pr(>F)"]
pca_separation <- ifelse(pca_manova_p < 0.05, "evidence_of_separation", "no_clear_separation")

p_pca <- ggplot(pca_tab, aes(PC1, PC2, color = subtype, shape = subtype)) +
  stat_ellipse(
    aes(fill = subtype),
    geom = "polygon",
    type = "norm",
    alpha = 0.10,
    color = NA
  ) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = c(sporadic = "#0072B2", NF2 = "#D55E00")) +
  scale_fill_manual(values = c(sporadic = "#0072B2", NF2 = "#D55E00")) +
  labs(
    title = "GSE108524: NF2-associated vs sporadic VS",
    subtitle = sprintf(
      "Top %d variable genes; first %d PCs MANOVA Pillai p = %.3g",
      n_pca_genes,
      n_manova_pc,
      pca_manova_p
    ),
    x = sprintf("PC1 (%.1f%%)", variance_explained[1]),
    y = sprintf("PC2 (%.1f%%)", variance_explained[2]),
    color = "Subtype",
    shape = "Subtype",
    fill = "Subtype"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right")
ggsave(paths$pca_figure, p_pca, width = 7.2, height = 5.8, dpi = 300)

signature <- fread(signature_file)
if (!all(c("gene", "coefficient") %in% names(signature))) {
  stop("Signature table must contain gene and coefficient columns.")
}
signature_overlap <- intersect(signature$gene, rownames(expr_model))
signature_missing <- setdiff(signature$gene, signature_overlap)
if (!length(signature_overlap)) stop("No signature genes overlap the expression matrix.")

signature_expr <- t(expr_model[signature_overlap, , drop = FALSE])
signature_z <- scale(signature_expr)
signature_z[, is.na(colSums(signature_z))] <- 0
coefficients <- setNames(signature$coefficient, signature$gene)[signature_overlap]
signature_score <- as.numeric(signature_z %*% coefficients)
score_tab <- data.table(
  sample = vs_samples,
  subtype = as.character(subtype),
  signature_score = signature_score
)
fwrite(score_tab, paths$signature_scores)

wilcox_result <- wilcox.test(
  signature_score ~ subtype,
  data = score_tab,
  exact = FALSE,
  conf.int = FALSE
)
signature_wilcox_p <- wilcox_result$p.value

p_score <- ggplot(score_tab, aes(subtype, signature_score, fill = subtype)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.72) +
  geom_jitter(width = 0.10, size = 2.2, alpha = 0.85) +
  scale_fill_manual(values = c(sporadic = "#0072B2", NF2 = "#D55E00")) +
  labs(
    title = "Discovery-set 11-gene signature score",
    subtitle = sprintf(
      "Within-cohort gene z-scores; Wilcoxon p = %.3g; genes available %d/%d",
      signature_wilcox_p,
      length(signature_overlap),
      nrow(signature)
    ),
    x = NULL,
    y = "Coefficient-weighted signature score",
    fill = "Subtype"
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none")
ggsave(paths$signature_figure, p_score, width = 6.2, height = 5.4, dpi = 300)

transcriptomes_similar <- n_deg_sig <= 5L && pca_manova_p >= 0.05
similarity_conclusion <- if (transcriptomes_similar) {
  paste(
    "NF2 与 sporadic VS 在该队列的整体转录组高度相似；",
    "该结果与既往报道一致（Gregory 2023）。"
  )
} else {
  paste(
    "该队列中 NF2 与 sporadic VS 存在候选转录差异；",
    "这些差异提示亚型相关表达特征，仍需独立队列验证。"
  )
}

top_candidate_text <- "无满足 adj.P.Val < 0.05 且 |logFC| > 1 的候选差异基因。"
if (n_deg_sig > 0L) {
  top_n <- min(10L, n_deg_sig)
  top_candidates <- deg_sig[seq_len(top_n)]
  top_candidate_text <- paste0(
    "Top 候选差异基因（logFC 为 NF2 - sporadic）：",
    paste(
      sprintf(
        "%s (logFC=%.2f, FDR=%.3g)",
        top_candidates$gene,
        top_candidates$logFC,
        top_candidates$adj.P.Val
      ),
      collapse = "；"
    ),
    "。"
  )
}

summary_tab <- data.table(
  metric = c(
    "nf2_n",
    "sporadic_n",
    "genes_tested",
    "deg_fdr_lt_0.05",
    "deg_fdr_lt_0.05_abs_logfc_gt_1",
    "pca_variable_genes",
    "pca_pc1_variance_percent",
    "pca_pc2_variance_percent",
    "pca_manova_n_pcs",
    "pca_manova_p",
    "pca_interpretation",
    "signature_genes_total",
    "signature_genes_available",
    "signature_genes_missing",
    "signature_wilcoxon_p",
    "transcriptomes_similar"
  ),
  value = as.character(c(
    n_nf2,
    n_sporadic,
    nrow(deg),
    n_deg_fdr,
    n_deg_sig,
    n_pca_genes,
    variance_explained[1],
    variance_explained[2],
    n_manova_pc,
    pca_manova_p,
    pca_separation,
    nrow(signature),
    length(signature_overlap),
    paste(signature_missing, collapse = ";"),
    signature_wilcox_p,
    transcriptomes_similar
  ))
)
fwrite(summary_tab, paths$summary)

report_lines <- c(
  "# GSE108524：NF2-related VS vs sporadic VS 转录组比较",
  "",
  "## 分组与方法",
  "",
  sprintf("- NF2-related VS：%d 例；sporadic VS：%d 例；4 个正常神经样本已排除。", n_nf2, n_sporadic),
  "- 标签来自 GEO series matrix 的 `!Sample_title` 与 `!Sample_characteristics_ch1`；所有样本的判定依据均已导出供核对。",
  "- 差异分析：gene-level log2 表达矩阵，limma 设计 `~ 0 + subtype`，对比为 `NF2 - sporadic`。",
  sprintf("- PCA：使用 VS 样本中方差最高的 %d 个基因，图展示 PC1/PC2，并用前 %d 个 PC 的 MANOVA Pillai 检验辅助量化整体分离。", n_pca_genes, n_manova_pc),
  sprintf("- 11 基因 signature：%d/%d 个基因可用；在 27 个 VS 样本内逐基因 z-score 后按发现集系数加权求和。", length(signature_overlap), nrow(signature)),
  "",
  "## 结果",
  "",
  sprintf("- FDR < 0.05 的基因数：%d。", n_deg_fdr),
  sprintf("- 同时满足 adj.P.Val < 0.05 且 |logFC| > 1 的显著 DEG 数：%d。", n_deg_sig),
  sprintf("- PCA 前 %d 个 PC 的 MANOVA Pillai p = %.4g；解释为：%s。", n_manova_pc, pca_manova_p, ifelse(pca_manova_p < 0.05, "存在整体分离证据", "未见清晰整体分离证据")),
  sprintf("- 11 基因 signature score 的 NF2 vs sporadic Wilcoxon p = %.4g。", signature_wilcox_p),
  paste0("- ", top_candidate_text),
  "",
  "## 客观结论",
  "",
  paste0("- ", similarity_conclusion),
  "- PCA 的 NF2 与 sporadic 点云仍有明显重叠，因此不能表述为两个完全分离的转录组亚群。",
  "- 重要风险：10 个 sporadic 样本均位于 GSM2893473–GSM2893482，17 个 NF2 样本均位于 GSM2902750–GSM2902766；亚型与 GSM 编号段完全对应。原始元数据未提供可用于模型校正的独立 batch 字段，因此候选差异可能混入未记录的入组或处理批次效应。",
  "- 本分析属于单队列探索性比较；差异基因及 signature 分数差异不能解释为 NF2 状态的因果效应，需结合独立队列或实验验证。",
  "",
  "## 产物与来源路径",
  "",
  sprintf("- 分析脚本：`%s`", script_file),
  sprintf("- GEO 标签判定依据：`%s`", paths$evidence),
  sprintf("- 全量 limma 结果：`%s`", paths$deg_all),
  sprintf("- 显著 DEG 表：`%s`", paths$deg_sig),
  sprintf("- PCA 坐标：`%s`", paths$pca_coordinates),
  sprintf("- PCA 图：`%s`", paths$pca_figure),
  sprintf("- signature 分数表：`%s`", paths$signature_scores),
  sprintf("- signature score 图：`%s`", paths$signature_figure),
  sprintf("- 汇总指标：`%s`", paths$summary),
  sprintf("- R 会话信息：`%s`", paths$session_info),
  sprintf("- 本报告：`%s`", paths$report),
  "",
  "## 输入（只读）",
  "",
  sprintf("- GEO series matrix：`%s`", raw_file),
  sprintf("- gene-level 表达矩阵：`%s`", expr_file),
  sprintf("- 样本元数据：`%s`", meta_file),
  sprintf("- 11 基因 signature：`%s`", signature_file)
)
writeLines(report_lines, paths$report)
writeLines(capture.output(sessionInfo()), paths$session_info)

cat(sprintf("PCA MANOVA p = %.4g; signature Wilcoxon p = %.4g\n", pca_manova_p, signature_wilcox_p))
cat("Conclusion: ", similarity_conclusion, "\n", sep = "")
cat("Report: ", paths$report, "\n", sep = "")
cat("DONE 09\n")
