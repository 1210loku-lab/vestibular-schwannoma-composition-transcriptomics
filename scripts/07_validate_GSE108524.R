#!/usr/bin/env Rscript

root <- Sys.getenv("PROJ_ROOT", getwd())
project_lib <- file.path(root, "data/processed/R_library")
dir.create(project_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(project_lib, .libPaths()))

if (!requireNamespace("hta20transcriptcluster.db", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install(
    "hta20transcriptcluster.db",
    lib = project_lib,
    update = FALSE,
    ask = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(AnnotationDbi)
  library(hta20transcriptcluster.db)
  library(ggplot2)
  library(ggrepel)
  library(pROC)
})

options(stringsAsFactors = FALSE)
set.seed(20260707)

raw_file <- file.path(root, "data/raw/GSE108524_matrix.txt.gz")
processed_dir <- file.path(root, "data/processed")
outdir <- file.path(root, "results/validation")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

lines <- readLines(gzfile(raw_file))
get_series_row <- function(prefix) {
  hit <- lines[startsWith(lines, prefix)]
  if (!length(hit)) stop(sprintf("Missing series-matrix row: %s", prefix))
  gsub('"', "", strsplit(hit[1], "\t", fixed = TRUE)[[1]][-1])
}

sample_ids <- get_series_row("!Sample_geo_accession")
characteristic_lines <- lines[startsWith(lines, "!Sample_characteristics_ch1")]
characteristics <- lapply(
  characteristic_lines,
  function(z) gsub('"', "", strsplit(z, "\t", fixed = TRUE)[[1]][-1])
)
disease_idx <- which(vapply(
  characteristics,
  function(z) any(grepl("^disease state:", z, ignore.case = TRUE)),
  logical(1)
))[1]
if (is.na(disease_idx)) stop("No disease-state characteristic found.")
disease_state <- sub(
  "^disease state:\\s*",
  "",
  characteristics[[disease_idx]],
  ignore.case = TRUE
)
group <- ifelse(
  grepl("^control$", disease_state, ignore.case = TRUE),
  "Control",
  "Tumor"
)

table_begin <- grep("^!series_matrix_table_begin", lines)
table_end <- grep("^!series_matrix_table_end", lines)
stopifnot(length(table_begin) == 1, length(table_end) == 1, table_end > table_begin)
expr_probe_dt <- fread(
  text = paste(lines[(table_begin + 1):(table_end - 1)], collapse = "\n"),
  header = TRUE
)
setnames(expr_probe_dt, 1, "probe")
expr_probe <- as.matrix(expr_probe_dt[, -1])
rownames(expr_probe) <- expr_probe_dt$probe
colnames(expr_probe) <- gsub('"', "", colnames(expr_probe))
storage.mode(expr_probe) <- "double"

stopifnot(
  identical(colnames(expr_probe), sample_ids),
  length(group) == ncol(expr_probe),
  sum(group == "Control") == 4,
  sum(group == "Tumor") == 27
)

was_log2_transformed <- FALSE
if (max(expr_probe, na.rm = TRUE) > 50) {
  expr_probe <- log2(expr_probe + 1)
  was_log2_transformed <- TRUE
}

probe_map <- AnnotationDbi::select(
  hta20transcriptcluster.db,
  keys = rownames(expr_probe),
  columns = "SYMBOL",
  keytype = "PROBEID"
)
probe_map <- unique(probe_map[
  !is.na(probe_map$SYMBOL) & nzchar(probe_map$SYMBOL),
  c("PROBEID", "SYMBOL")
])
probe_map$mean_expression <- rowMeans(
  expr_probe[probe_map$PROBEID, , drop = FALSE],
  na.rm = TRUE
)
probe_map <- probe_map[
  order(probe_map$SYMBOL, -probe_map$mean_expression, probe_map$PROBEID),
]
probe_map <- probe_map[!duplicated(probe_map$SYMBOL), ]

expr_gene <- expr_probe[probe_map$PROBEID, , drop = FALSE]
rownames(expr_gene) <- probe_map$SYMBOL
expr_gene <- expr_gene[!duplicated(rownames(expr_gene)), , drop = FALSE]

validation_meta <- data.table(
  sample = sample_ids,
  group = group,
  disease_state = disease_state
)
fwrite(
  data.table(gene = rownames(expr_gene), expr_gene),
  file.path(processed_dir, "GSE108524_expr_gene.csv")
)
fwrite(
  validation_meta,
  file.path(processed_dir, "GSE108524_meta.csv")
)

hub <- fread(file.path(root, "results/hub/hub_all.csv"))
hub_overlap <- intersect(hub$gene, rownames(expr_gene))
validation_logfc <- rowMeans(
  expr_gene[, group == "Tumor", drop = FALSE],
  na.rm = TRUE
) - rowMeans(
  expr_gene[, group == "Control", drop = FALSE],
  na.rm = TRUE
)

hub_concordance <- merge(
  hub[, .(
    gene,
    module,
    direction,
    discovery_logFC = logFC
  )],
  data.table(
    gene = names(validation_logfc),
    validation_logFC = as.numeric(validation_logfc)
  ),
  by = "gene",
  all.x = TRUE
)
hub_concordance[, direction_concordant := (
  sign(discovery_logFC) == sign(validation_logFC)
)]
fwrite(hub_concordance, file.path(outdir, "hub_logFC_concordance.csv"))

hub_complete <- hub_concordance[
  !is.na(discovery_logFC) & !is.na(validation_logFC)
]
stopifnot(nrow(hub_complete) >= 3)
concordant_n <- sum(hub_complete$direction_concordant)
concordance_rate <- mean(hub_complete$direction_concordant)
concordance_test <- binom.test(
  concordant_n,
  nrow(hub_complete),
  p = 0.5,
  alternative = "greater"
)
pearson_test <- cor.test(
  hub_complete$discovery_logFC,
  hub_complete$validation_logFC,
  method = "pearson"
)
spearman_test <- suppressWarnings(cor.test(
  hub_complete$discovery_logFC,
  hub_complete$validation_logFC,
  method = "spearman",
  exact = FALSE
))

label_data <- hub_complete[
  order(-abs(discovery_logFC * validation_logFC))
][seq_len(min(12, nrow(hub_complete)))]
p_concordance <- ggplot(
  hub_complete,
  aes(x = discovery_logFC, y = validation_logFC, color = module)
) +
  geom_hline(yintercept = 0, color = "grey65", linewidth = 0.4) +
  geom_vline(xintercept = 0, color = "grey65", linewidth = 0.4) +
  geom_point(alpha = 0.75, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.7) +
  geom_text_repel(
    data = label_data,
    aes(label = gene),
    size = 3,
    max.overlaps = 20,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c(turquoise = "#009E73", yellow = "#E69F00")) +
  labs(
    title = "Hub-gene logFC concordance across datasets",
    subtitle = sprintf(
      "Overlap %d/%d; direction concordance %.1f%%; Pearson r = %.3f",
      nrow(hub_complete),
      nrow(hub),
      100 * concordance_rate,
      unname(pearson_test$estimate)
    ),
    x = "GSE39645 log2FC (VS vs control)",
    y = "GSE108524 log2FC (VS vs control)",
    color = "Module"
  ) +
  theme_bw(base_size = 11)
ggsave(
  file.path(outdir, "fig_logFC_concordance.png"),
  p_concordance,
  width = 7,
  height = 5.8,
  dpi = 180
)

signature <- fread(file.path(root, "results/signature/signature_genes.csv"))
model_bundle <- readRDS(file.path(root, "results/signature/model.rds"))
signature_overlap <- intersect(signature$gene, rownames(expr_gene))
signature_missing <- setdiff(signature$gene, rownames(expr_gene))

fwrite(
  data.table(
    gene = signature$gene,
    present_in_GSE108524 = signature$gene %in% signature_overlap
  ),
  file.path(outdir, "signature_gene_overlap.csv")
)
fwrite(
  data.table(gene = signature_missing),
  file.path(outdir, "signature_missing_genes.csv")
)
stopifnot(length(signature_overlap) >= 1)

validation_sig <- t(expr_gene[signature_overlap, , drop = FALSE])
validation_sig_z <- scale(validation_sig)
validation_sig_z[, is.na(colSums(validation_sig_z))] <- 0
coef_overlap <- model_bundle$coefficients[signature_overlap]
signature_score <- as.numeric(validation_sig_z %*% coef_overlap)

score_tab <- data.table(
  sample = sample_ids,
  group = group,
  disease_state = disease_state,
  signature_score = signature_score
)
fwrite(score_tab, file.path(outdir, "signature_scores.csv"))

roc_validation <- suppressWarnings(roc(
  response = factor(group, levels = c("Control", "Tumor")),
  predictor = signature_score,
  levels = c("Control", "Tumor"),
  direction = "<",
  quiet = TRUE,
  ci = TRUE
))
validation_auc <- as.numeric(auc(roc_validation))
validation_auc_ci <- as.numeric(suppressWarnings(ci.auc(roc_validation)))

png(
  file.path(outdir, "fig_roc_validation.png"),
  width = 1000,
  height = 850,
  res = 150
)
plot(
  roc_validation,
  col = "#D55E00",
  lwd = 3,
  legacy.axes = TRUE,
  identity = TRUE,
  identity.col = "grey60",
  identity.lty = 2,
  main = "GSE108524 external validation ROC"
)
legend(
  "bottomright",
  legend = sprintf(
    "AUC = %.3f (95%% CI %.3f\u2013%.3f)\nGenes available: %d/%d",
    validation_auc,
    validation_auc_ci[1],
    validation_auc_ci[3],
    length(signature_overlap),
    nrow(signature)
  ),
  bty = "n"
)
dev.off()

validation_summary <- data.table(
  metric = c(
    "validation_samples_total",
    "validation_controls",
    "validation_tumors",
    "validation_genes_after_annotation",
    "log2_transformation_applied",
    "hub_genes_total",
    "hub_genes_validation_overlap",
    "logFC_direction_concordance",
    "logFC_concordance_binomial_p",
    "logFC_pearson_r",
    "logFC_pearson_p",
    "logFC_spearman_rho",
    "logFC_spearman_p",
    "signature_genes_total",
    "signature_genes_validation_overlap",
    "signature_genes_missing",
    "signature_validation_auc",
    "signature_validation_auc_ci_low",
    "signature_validation_auc_ci_high"
  ),
  value = as.character(c(
    ncol(expr_gene),
    sum(group == "Control"),
    sum(group == "Tumor"),
    nrow(expr_gene),
    was_log2_transformed,
    nrow(hub),
    nrow(hub_complete),
    concordance_rate,
    concordance_test$p.value,
    unname(pearson_test$estimate),
    pearson_test$p.value,
    unname(spearman_test$estimate),
    spearman_test$p.value,
    nrow(signature),
    length(signature_overlap),
    length(signature_missing),
    validation_auc,
    validation_auc_ci[1],
    validation_auc_ci[3]
  ))
)
fwrite(validation_summary, file.path(outdir, "validation_summary.csv"))

cat(sprintf(
  paste0(
    "GSE108524: %d controls, %d tumors, %d annotated genes\n",
    "Hub overlap: %d/%d; direction concordance: %.2f%%; ",
    "Pearson r = %.4f; Spearman rho = %.4f\n",
    "Signature overlap: %d/%d; missing: %s\n",
    "External validation AUC: %.4f (95%% CI %.4f-%.4f)\n"
  ),
  sum(group == "Control"),
  sum(group == "Tumor"),
  nrow(expr_gene),
  nrow(hub_complete),
  nrow(hub),
  100 * concordance_rate,
  unname(pearson_test$estimate),
  unname(spearman_test$estimate),
  length(signature_overlap),
  nrow(signature),
  ifelse(length(signature_missing), paste(signature_missing, collapse = ", "), "none"),
  validation_auc,
  validation_auc_ci[1],
  validation_auc_ci[3]
))
cat("DONE 07\n")
