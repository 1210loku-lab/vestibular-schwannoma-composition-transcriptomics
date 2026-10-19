#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

ROOT <- Sys.getenv("VS_INPUT_ROOT", unset = ".")
EXPORT_ROOT <- Sys.getenv("VS_OUTPUT_ROOT", unset = ".")
ANALYSIS <- file.path(EXPORT_ROOT, "outputs")
FIGURES <- file.path(EXPORT_ROOT, "figures")
SCRIPT_DIR <- file.path(EXPORT_ROOT, "scripts")
MANUSCRIPT <- file.path(EXPORT_ROOT, "notes")
COUNTS_PATH <- file.path(ROOT, "data", "processed", "GSE230375_pseudobulk_counts.csv")
CELL_META_PATH <- file.path(ROOT, "data", "processed", "GSE230375_cell_meta.csv")
DONOR_META_PATH <- file.path(ROOT, "results", "pseudobulk", "pseudobulk_sample_celltype_ncells.csv")

dir.create(ANALYSIS, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURES, recursive = TRUE, showWarnings = FALSE)
dir.create(SCRIPT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(MANUSCRIPT, recursive = TRUE, showWarnings = FALSE)

log_note <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), paste0(...)))
  flush.console()
}

required_packages <- c("DESeq2", "apeglm", "ggplot2", "patchwork", "svglite", "ragg")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

lineage_map <- data.frame(
  cell_type = c("Schwann/肿瘤细胞", "巨噬/髓系", "T细胞", "成纤维", "内皮", "周细胞"),
  lineage = c("Schwann/tumor-associated annotation (combined)", "Myeloid", "T cell", "Fibroblast", "Endothelial", "Pericyte"),
  lineage_id = c("Schwann_tumor_associated", "Myeloid", "T_cell", "Fibroblast", "Endothelial", "Pericyte"),
  stringsAsFactors = FALSE
)
donor_order <- c("N1", "N2", paste0("T", 1:7))
fixed_genes <- c("NRG1", "L1CAM", "FCGBP", "PRRX1")
selected_panels <- data.frame(
  gene = fixed_genes,
  cell_type = c("Schwann/肿瘤细胞", "Schwann/肿瘤细胞", "巨噬/髓系", "成纤维"),
  stringsAsFactors = FALSE
)

log_note("donor-reference starting; reading processed pseudobulk counts and cell metadata. The Seurat RDS is not read.")
counts_frame <- read.csv(COUNTS_PATH, check.names = FALSE, stringsAsFactors = FALSE)
genes <- as.character(counts_frame[[1]])
count_mat <- as.matrix(counts_frame[, -1, drop = FALSE])
storage.mode(count_mat) <- "double"
rownames(count_mat) <- genes
cell_meta <- read.csv(CELL_META_PATH, check.names = FALSE, stringsAsFactors = FALSE)
donor_meta <- read.csv(DONOR_META_PATH, check.names = FALSE, stringsAsFactors = FALSE)

if (anyNA(genes) || any(!nzchar(genes))) stop("Gene identifiers contain missing/empty values")
if (anyDuplicated(genes)) stop("Gene identifiers are not unique")
if (any(!is.finite(count_mat)) || any(count_mat < 0)) stop("Pseudobulk counts must be finite and nonnegative")
if (any(abs(count_mat - round(count_mat)) > 1e-8)) stop("Pseudobulk counts are not integer-valued")
if (anyDuplicated(colnames(count_mat))) stop("Pseudobulk column identifiers are not unique")
if (!all(c("barcode", "sample", "group", "cell_type") %in% names(cell_meta))) {
  stop("Cell metadata is missing required fields")
}
if (!all(cell_meta$cell_type %in% lineage_map$cell_type)) {
  stop("Unexpected cell type labels in processed metadata")
}

donors <- unique(cell_meta[, c("sample", "group"), drop = FALSE])
if (nrow(donors) != length(donor_order) || anyDuplicated(donors$sample)) {
  stop("Expected one group label for each of nine donors")
}
if (!setequal(as.character(donors$sample), donor_order)) stop("Unexpected donor IDs in cell metadata")
donors <- donors[match(donor_order, donors$sample), , drop = FALSE]
if (!identical(as.character(donors$group), c("Nerve", "Nerve", rep("VS", 7)))) {
  stop("Donor group labels do not match the expected two nerve and seven VS donors")
}

meta_counts <- aggregate(
  x = list(n_cells = rep.int(1L, nrow(cell_meta))),
  by = list(sample = as.character(cell_meta$sample), cell_type = as.character(cell_meta$cell_type)),
  FUN = sum
)
meta_counts$id <- paste(meta_counts$sample, meta_counts$cell_type, sep = "@@")
expected_ids <- as.character(meta_counts$id)
if (length(expected_ids) != 54L || anyDuplicated(expected_ids)) stop("Metadata does not define 54 unique donor-compartment pairs")
if (!setequal(colnames(count_mat), expected_ids)) stop("Pseudobulk matrix IDs do not exactly match metadata-derived donor-compartment IDs")
if (nrow(donor_meta) != 54L || anyDuplicated(donor_meta$id)) stop("Existing donor metadata is not 54 unique pairs")

metadata_check <- merge(
  meta_counts[, c("id", "sample", "cell_type", "n_cells")],
  donor_meta[, c("id", "sample", "cell_type", "n_cells", "group")],
  by = "id", suffixes = c("_from_cell_metadata", "_existing"), all = TRUE, sort = FALSE
)
if (nrow(metadata_check) != 54L || anyNA(metadata_check)) stop("Existing donor metadata does not align with cell metadata")
if (any(metadata_check$n_cells_from_cell_metadata != metadata_check$n_cells_existing) ||
    any(metadata_check$sample_from_cell_metadata != metadata_check$sample_existing) ||
    any(metadata_check$cell_type_from_cell_metadata != metadata_check$cell_type_existing)) {
  stop("Independent cell metadata recount differs from the existing donor metadata")
}
sample_group <- setNames(as.character(donors$group), as.character(donors$sample))
metadata_check$group_from_cell_metadata <- unname(sample_group[metadata_check$sample_from_cell_metadata])
if (any(metadata_check$group_from_cell_metadata != metadata_check$group_existing)) {
  stop("Existing donor group labels differ from cell metadata")
}
if (any(metadata_check$n_cells_from_cell_metadata < 10L)) stop("A donor-compartment is below the original 10-cell inclusion threshold")

lineage_for <- setNames(lineage_map$lineage, lineage_map$cell_type)
lineage_id_for <- setNames(lineage_map$lineage_id, lineage_map$cell_type)
group_label <- c(Nerve = "Nerve control", VS = "Vestibular schwannoma (VS)")

donor_table <- data.frame(
  donor_id = as.character(donors$sample),
  group = unname(group_label[as.character(donors$group)]),
  stringsAsFactors = FALSE
)
for (ct in lineage_map$cell_type) {
  ct_rows <- meta_counts[meta_counts$cell_type == ct, , drop = FALSE]
  donor_table[[unname(lineage_id_for[ct])]] <- ct_rows$n_cells[match(donor_table$donor_id, ct_rows$sample)]
}
write.csv(donor_table, file.path(ANALYSIS, "donor_donor_compartment_cell_counts.csv"), row.names = FALSE, na = "")

group_count_rows <- do.call(rbind, lapply(lineage_map$cell_type, function(ct) {
  do.call(rbind, lapply(c("Nerve", "VS"), function(g) {
    data.frame(
      cell_type = ct,
      lineage = unname(lineage_for[ct]),
      group = unname(group_label[g]),
      biological_donors = sum(donors$group == g),
      donors_meeting_n_cells_ge_10 = sum(meta_counts$n_cells[meta_counts$cell_type == ct &
        meta_counts$sample %in% donors$sample[donors$group == g]] >= 10L),
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(group_count_rows, file.path(ANALYSIS, "donor_compartment_group_donor_counts.csv"), row.names = FALSE)

matrix_meta <- data.frame(
  id = colnames(count_mat),
  sample = sub("@@.*$", "", colnames(count_mat)),
  cell_type = sub("^.*@@", "", colnames(count_mat)),
  stringsAsFactors = FALSE
)
matrix_meta$n_cells <- meta_counts$n_cells[match(matrix_meta$id, meta_counts$id)]
matrix_meta$group <- unname(sample_group[matrix_meta$sample])
matrix_meta$lineage <- unname(lineage_for[matrix_meta$cell_type])
matrix_meta$lineage_id <- unname(lineage_id_for[matrix_meta$cell_type])
matrix_meta$count_sum <- colSums(count_mat)
matrix_meta$nonzero_genes <- colSums(count_mat > 0)
matrix_meta$passed_original_n_cells_gate <- matrix_meta$n_cells >= 10L
write.csv(matrix_meta, file.path(ANALYSIS, "donor_pseudobulk_column_qc.csv"), row.names = FALSE)

if (length(fixed_genes) != 4L || !all(fixed_genes %in% genes)) stop("One or more predeclared representative genes are absent")
if (ncol(count_mat) != 54L || nrow(count_mat) != length(unique(genes))) stop("Unexpected count matrix dimensions")
log_note("Input checks passed: genes=", nrow(count_mat), "; donor-compartment columns=", ncol(count_mat),
         "; total pseudobulk counts=", format(sum(count_mat), scientific = FALSE),
         "; all 54 pairs pass n_cells>=10; independently confirmed Schwann nerve donors N1=",
         meta_counts$n_cells[meta_counts$sample == "N1" & meta_counts$cell_type == "Schwann/肿瘤细胞"],
         " and N2=", meta_counts$n_cells[meta_counts$sample == "N2" & meta_counts$cell_type == "Schwann/肿瘤细胞"])

all_effects <- vector("list", length(lineage_map$cell_type))
names(all_effects) <- lineage_map$cell_type
expression_rows <- vector("list", length(lineage_map$cell_type))
names(expression_rows) <- lineage_map$cell_type

for (ct in lineage_map$cell_type) {
  ids <- paste(donor_order, ct, sep = "@@")
  if (!all(ids %in% colnames(count_mat))) stop("Missing donor columns for compartment: ", ct)
  cd <- data.frame(
    group = factor(unname(sample_group[donor_order]), levels = c("Nerve", "VS")),
    row.names = ids
  )
  m <- count_mat[, ids, drop = FALSE]
  keep <- rowSums(m) >= 10
  if (sum(keep) == 0L) stop("No genes pass the original total-count filter for ", ct)
  dds <- DESeqDataSetFromMatrix(countData = round(m[keep, , drop = FALSE]), colData = cd, design = ~ group)
  dds <- DESeq(dds, quiet = TRUE)
  coef_name <- "group_VS_vs_Nerve"
  if (!coef_name %in% resultsNames(dds)) {
    stop("Expected DESeq2 coefficient is absent for ", ct, ": ", paste(resultsNames(dds), collapse = ","))
  }
  raw <- results(dds, contrast = c("group", "VS", "Nerve"), independentFiltering = TRUE)
  shrink_result <- lfcShrink(dds, coef = coef_name, type = "apeglm", returnList = TRUE)
  shr_df <- as.data.frame(shrink_result$res)
  raw_df <- as.data.frame(raw)
  if (!identical(rownames(raw_df), rownames(shr_df))) {
    shr_df <- shr_df[match(rownames(raw_df), rownames(shr_df)), , drop = FALSE]
  }
  raw_beta_conv <- mcols(dds)$betaConv
  if (is.null(raw_beta_conv) || length(raw_beta_conv) != nrow(dds)) {
    stop("DESeq2 betaConv diagnostic is missing or has the wrong length for ", ct)
  }
  apeglm_diag <- shrink_result$fit$diag
  if (is.null(dim(apeglm_diag)) || nrow(apeglm_diag) != nrow(dds) ||
      !"conv" %in% colnames(apeglm_diag)) {
    stop("apeglm convergence diagnostic is missing or malformed for ", ct)
  }
  raw_beta_conv_by_gene <- setNames(as.logical(raw_beta_conv), rownames(dds))
  apeglm_conv_by_gene <- setNames(as.numeric(apeglm_diag[, "conv"]), rownames(dds))
  effects <- data.frame(
    gene = rownames(raw_df),
    cell_type = ct,
    lineage = unname(lineage_for[ct]),
    n_VS = sum(cd$group == "VS"),
    n_Nerve = sum(cd$group == "Nerve"),
    total_raw_count = rowSums(m[keep, , drop = FALSE])[rownames(raw_df)],
    baseMean = raw_df$baseMean,
    raw_log2FC = raw_df$log2FoldChange,
    raw_WaldSE = raw_df$lfcSE,
    raw_pvalue = raw_df$pvalue,
    raw_padj = raw_df$padj,
    raw_beta_converged = unname(raw_beta_conv_by_gene[rownames(raw_df)]),
    apeglm_shrunken_log2FC = shr_df$log2FoldChange,
    apeglm_posterior_SD = shr_df$lfcSE,
    apeglm_convergence_code = unname(apeglm_conv_by_gene[rownames(raw_df)]),
    apeglm_converged = unname(apeglm_conv_by_gene[rownames(raw_df)] == 0),
    stringsAsFactors = FALSE
  )
  all_effects[[ct]] <- effects
  norm <- sweep(m[fixed_genes, ids, drop = FALSE], 2, sizeFactors(dds), "/")
  expression_rows[[ct]] <- do.call(rbind, lapply(fixed_genes, function(gene) {
    cell_counts <- meta_counts$n_cells[match(ids, meta_counts$id)]
    data.frame(
      gene = gene,
      cell_type = ct,
      lineage = unname(lineage_for[ct]),
      donor_id = donor_order,
      group = unname(group_label[unname(sample_group[donor_order])]),
      group_code = unname(sample_group[donor_order]),
      n_cells = cell_counts,
      DESeq2_normalized_count = as.numeric(norm[gene, ids]),
      log2_normalized_count_plus_1 = log2(as.numeric(norm[gene, ids]) + 1),
      stringsAsFactors = FALSE
    )
  }))
  log_note("DESeq2 fit complete: ", unname(lineage_for[ct]), "; tested genes=", nrow(effects),
           "; Nerve donors=", sum(cd$group == "Nerve"), "; VS donors=", sum(cd$group == "VS"),
           "; apeglm=", packageVersion("apeglm"))
}

effects_long <- do.call(rbind, all_effects)
rownames(effects_long) <- NULL
effects_path <- file.path(ANALYSIS, "donor_deseq2_gene_effects_all_compartments.csv.gz")
write.csv(effects_long, gzfile(effects_path, "wt"), row.names = FALSE, na = "")

fit_convergence_summary <- do.call(rbind, lapply(lineage_map$cell_type, function(ct) {
  x <- effects_long[effects_long$cell_type == ct, , drop = FALSE]
  data.frame(
    cell_type = ct,
    lineage = unname(lineage_for[ct]),
    genes_tested = nrow(x),
    raw_beta_converged = sum(x$raw_beta_converged %in% TRUE),
    raw_beta_not_converged = sum(x$raw_beta_converged %in% FALSE),
    raw_beta_status_missing = sum(is.na(x$raw_beta_converged)),
    apeglm_converged = sum(x$apeglm_converged %in% TRUE),
    apeglm_not_converged = sum(x$apeglm_converged %in% FALSE),
    apeglm_convergence_status_missing = sum(is.na(x$apeglm_converged)),
    apeglm_posterior_SD_missing = sum(is.na(x$apeglm_posterior_SD)),
    stringsAsFactors = FALSE
  )
}))
write.csv(fit_convergence_summary, file.path(ANALYSIS, "donor_fit_convergence_summary.csv"), row.names = FALSE)

expression_table <- do.call(rbind, expression_rows)
rownames(expression_table) <- NULL
expression_table$gene <- factor(expression_table$gene, levels = fixed_genes)
expression_table$lineage <- factor(expression_table$lineage, levels = lineage_map$lineage)
expression_table$group <- factor(expression_table$group,
  levels = c("Nerve control", "Vestibular schwannoma (VS)"))
expression_table$group_code <- factor(expression_table$group_code, levels = c("Nerve", "VS"))
expression_table$cell_type <- as.character(expression_table$cell_type)
write.csv(expression_table, file.path(ANALYSIS, "donor_candidate_expression_donor_level.csv"), row.names = FALSE, na = "")

expression_summary <- do.call(rbind, lapply(split(expression_table,
  interaction(expression_table$gene, expression_table$cell_type, expression_table$group_code, drop = TRUE)), function(x) {
    data.frame(
      gene = as.character(x$gene[1]),
      cell_type = x$cell_type[1],
      lineage = as.character(x$lineage[1]),
      group = as.character(x$group[1]),
      n_donors = nrow(x),
      mean_log2_normalized_count_plus_1 = mean(x$log2_normalized_count_plus_1),
      SD_log2_normalized_count_plus_1 = if (nrow(x) > 1L) sd(x$log2_normalized_count_plus_1) else NA_real_,
      mean_DESeq2_normalized_count = mean(x$DESeq2_normalized_count),
      SD_DESeq2_normalized_count = if (nrow(x) > 1L) sd(x$DESeq2_normalized_count) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
rownames(expression_summary) <- NULL
write.csv(expression_summary, file.path(ANALYSIS, "donor_candidate_expression_summary.csv"), row.names = FALSE, na = "")

candidate_effects <- effects_long[effects_long$gene %in% fixed_genes, , drop = FALSE]
candidate_grid <- expand.grid(
  gene = fixed_genes,
  cell_type = lineage_map$cell_type,
  stringsAsFactors = FALSE
)
candidate_grid$lineage <- unname(lineage_for[candidate_grid$cell_type])
candidate_effects <- merge(candidate_grid, candidate_effects,
  by = c("gene", "cell_type", "lineage"), all.x = TRUE, sort = FALSE)
candidate_effects$tested_in_DESeq2 <- !is.na(candidate_effects$baseMean)
candidate_effects$fit_status <- ifelse(candidate_effects$tested_in_DESeq2,
  "tested_after_total_count_ge_10_filter", "not_tested_total_count_below_10")
candidate_effects$n_VS[is.na(candidate_effects$n_VS)] <- 7L
candidate_effects$n_Nerve[is.na(candidate_effects$n_Nerve)] <- 2L
candidate_effects$gene <- factor(candidate_effects$gene, levels = fixed_genes)
candidate_effects$cell_type <- factor(candidate_effects$cell_type, levels = lineage_map$cell_type)
candidate_effects$lineage <- factor(candidate_effects$lineage, levels = lineage_map$lineage)
candidate_effects <- candidate_effects[order(candidate_effects$gene, candidate_effects$cell_type), , drop = FALSE]
if (nrow(candidate_effects) != length(fixed_genes) * nrow(lineage_map)) stop("Candidate effect table does not contain all 24 gene-compartment pairs")
write.csv(candidate_effects, file.path(ANALYSIS, "donor_candidate_effects_all6_compartments.csv"), row.names = FALSE, na = "")
candidate_convergence <- candidate_effects[, c(
  "gene", "cell_type", "lineage", "fit_status", "n_VS", "n_Nerve",
  "raw_beta_converged", "raw_log2FC", "raw_WaldSE",
  "apeglm_converged", "apeglm_convergence_code",
  "apeglm_shrunken_log2FC", "apeglm_posterior_SD"
), drop = FALSE]
write.csv(candidate_convergence, file.path(ANALYSIS, "donor_candidate_convergence.csv"), row.names = FALSE, na = "")

plot_source <- expression_table[expression_table$gene %in% fixed_genes, , drop = FALSE]
selected_panels$lineage <- unname(lineage_for[selected_panels$cell_type])
panel_source <- do.call(rbind, lapply(seq_len(nrow(selected_panels)), function(i) {
  x <- plot_source[as.character(plot_source$gene) == selected_panels$gene[i] &
    plot_source$cell_type == selected_panels$cell_type[i], , drop = FALSE]
  x$panel_id <- letters[i]
  x$panel_gene <- selected_panels$gene[i]
  x$panel_lineage <- selected_panels$lineage[i]
  x
}))
write.csv(panel_source, file.path(ANALYSIS, "donor_Figure_S17_source.csv"), row.names = FALSE, na = "")

source(file.path(SCRIPT_DIR, "legacy_donor_expression_preview.R"), local = TRUE)

reference_table <- data.frame(
  reference_method = c("MuSiC", "NNLS"),
  GSE230375_input = c(
    "Raw-count SingleCellExperiment; colData carries sample and cell_type",
    "Processed cell-type mean log-normalized expression"
  ),
  construction = c(
    "MuSiC receives sample as the subject label and cell_type as the cluster label; cell-type profiles are informed by between-donor consistency weighting.",
    "Cross-cell mean expression per broad cell type from the processed matrix; pooled cell-level means, not equal-weight donor means."
  ),
  transformation_or_selection = c(
    "Raw counts are supplied to MuSiC; no marker-only restriction is applied in the original all-common-gene MuSiC analysis.",
    "Clamp mean log-normalized values at zero, apply expm1, rank cell-type specificity as log2((target+0.05)/(maximum other+0.05)), and union up to 200 markers per type with expression >0.05."
  ),
  boundary = c(
    "Shares the GSE230375 source with pseudobulk and is not independent validation.",
    "Marker-restricted NNLS provides a sensitivity reference; estimated fractions are not validated absolute tissue proportions."
  ),
  stringsAsFactors = FALSE
)
write.csv(reference_table, file.path(ANALYSIS, "donor_reference_construction.csv"), row.names = FALSE)

methods_text <- c(
  "### Donor-unit pseudobulk and reference construction",
  "",
  "For GSE230375, raw UMI counts were summed within each donor and broad cell compartment, retaining donor-compartment pseudobulks with at least 10 annotated cells. This yielded seven vestibular schwannoma donors and two nerve-control donors in each of six compartments: the broad Schwann/tumor-associated annotation, myeloid cells, T cells, fibroblasts, endothelial cells, and pericytes. The available annotation combines normal and tumor-associated Schwann states, so this compartment is not a pure neoplastic-cell population. The donor was the biological unit; cells were not treated as independent replicates. Per-donor captured-cell counts document provenance and threshold eligibility and are not estimates of tumor purity or weights used by MuSiC.",
  "",
  "Within each compartment, genes with summed counts below 10 were excluded, and DESeq2 was fit with a negative-binomial Wald model using group (VS versus nerve control) as the predictor. Size factors and normalized counts used the DESeq2 median-ratio procedure. The analysis reports the unshrunk log2 fold change, Wald standard error, Wald P value, and Benjamini–Hochberg adjusted P value separately from the apeglm-shrunken effect and its posterior standard deviation. Candidate expression displays use DESeq2-normalized counts as log2(count + 1); values are not rescaled to make the nerve-control mean equal to one.",
  "",
  "The GSE230375 MuSiC reference was constructed from a SingleCellExperiment containing raw counts and donor (`sample`) and cell-type (`cell_type`) labels, with the donor label supplied for MuSiC's between-donor weighting. The separate marker-restricted NNLS reference used pooled cross-cell means of log-normalized expression for each broad cell type, transformed with expm1 after truncation below zero; markers were selected by cell-type specificity, with up to 200 markers per type. This is not an equal-weight average of donor means. Normal and tumor-associated Schwann states remain merged. Because GSE230375 supplies both the pseudobulk analysis and the primary reference, these analyses are not independent validation. GSE216783 is a tumor-only second reference and cannot provide a VS-versus-nerve contrast. Agreement between references or recovery of expected ratios does not establish absolute tissue proportions or validate within-compartment disease effects."
)
writeLines(methods_text, file.path(MANUSCRIPT, "donor_reference_methods.md"), useBytes = TRUE)

log_note("Outputs written: donor table 9x6; group-count table 12 rows; column QC 54 rows; six DESeq2 fits; effects rows=",
         nrow(effects_long), "; expression donor rows=", nrow(expression_table),
         "; candidate effect pairs=", nrow(candidate_effects), "; raw beta convergence failures=",
         sum(effects_long$raw_beta_converged %in% FALSE), "; apeglm convergence failures=",
         sum(effects_long$apeglm_converged %in% FALSE), "; S17 PDF/PNG/SVG")
log_note("donor-reference analysis completed successfully.")
