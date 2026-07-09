#!/usr/bin/env Rscript
# 14_gse216783_second_ref.R
# 加深3 (TASK-D): GSE216783 (Barrett 2024, 15 sporadic VS scRNA/snRNA) 作为
#   (1) VS 内部独立的细胞区室定位 (candidate localization 复现);
#   (2) 第二个 tumour-derived 参考, 对 GSE39645/GSE108524 做反卷积敏感性分析。
# 不做 VS-vs-nerve pseudobulk; 不宣称破除单参考循环; 措辞: 候选/与…一致/提示, 不写因果。
# 重要: GSE216783 无正常神经对照, 且与已引的 Barrett 2024 同源 → 仅作敏感性/定位, 非独立生物学验证。
# 数据: data/raw/GSE216783/ 为 Cell Ranger RAW 矩阵 (36601 × 6.79M barcodes), 需 cell-calling。
# 本机 16GB: 每样本读入后立即按 QC 过滤到真实细胞再合并, 避免 OOM。

options(stringsAsFactors = FALSE, warn = 1)
suppressMessages({
  library(Matrix)
  library(ggplot2)
})

root <- Sys.getenv("PROJ_ROOT", getwd())
rawdir <- file.path(root, "data", "raw", "GSE216783")
outdir <- file.path(root, "results", "gse216783")
procd  <- file.path(root, "data", "processed")
docdir <- file.path(root, "docs")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(...) message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))

# ---- cell-calling QC 阈值 (raw 10x) ----
MIN_COUNTS  <- 500    # nCount 下限 (cell-calling)
MIN_FEAT    <- 200
MAX_MT_PCT  <- 20

MARKERS <- list(
  "Schwann/tumour" = c("SOX10","S100B","PMP22","MPZ","PLP1","NGFR"),
  "Myeloid"        = c("CD68","CD163","LYZ","AIF1"),
  "T cell"         = c("CD3D","CD3E","CD2"),
  "Fibroblast"     = c("DCN","LUM","COL1A1","PDGFRA"),
  "Endothelial"    = c("PECAM1","VWF","CLDN5"),
  "Pericyte"       = c("RGS5","ACTA2")
)
CELLTYPES <- names(MARKERS)
CANDIDATES <- c("NRG1","L1CAM","NCAM2","FCGBP","PRRX1","AR","SHOX2","HRH1",
                "RRM2B","SLC16A7","PDGFRA")

# GSM <-> SCH <-> assay
samples <- read.csv(file.path(outdir, "sample_metadata.csv"))
samples$sample <- sub("_.*$", "", samples$title)   # SCH1 ...

# ---- 读入单样本 RAW 10x, 立即 cell-calling ----
read_qc_sample <- function(gsm, sch, assay) {
  pre <- paste0(substr(gsm, 1, nchar(gsm) - 3), "nnn")
  mtx  <- file.path(rawdir, paste0(gsm, "_", sch, "-matrix.mtx.gz"))
  bc   <- file.path(rawdir, paste0(gsm, "_", sch, "-barcodes.tsv.gz"))
  ft   <- file.path(rawdir, paste0(gsm, "_", sch, "-features.tsv.gz"))
  if (!all(file.exists(mtx, bc, ft))) { log_msg("MISSING files for ", sch); return(NULL) }
  m <- tryCatch(readMM(mtx), error = function(e) { log_msg("  SKIP ", sch, " (readMM failed: ", conditionMessage(e), ")"); NULL })
  if (is.null(m)) return(NULL)
  m <- as(m, "CsparseMatrix")
  barcodes <- readLines(bc)
  feats <- read.delim(ft, header = FALSE)
  sym <- feats[[2]]
  # 基因符号去重 (保留首次)
  rownames(m) <- make.unique(sym)
  colnames(m) <- paste0(sch, "_", barcodes)
  # cell-calling: 列过滤
  cs <- Matrix::colSums(m)
  nf <- Matrix::colSums(m > 0)
  mt_idx <- grep("^MT-", rownames(m))
  mtpct <- if (length(mt_idx)) 100 * Matrix::colSums(m[mt_idx, , drop = FALSE]) / pmax(cs, 1) else rep(0, ncol(m))
  keep <- cs >= MIN_COUNTS & nf >= MIN_FEAT & mtpct < MAX_MT_PCT
  m <- m[, keep, drop = FALSE]
  log_msg("  ", sch, " [", assay, "] cells kept=", ncol(m), " / raw barcodes=", length(barcodes))
  list(mat = m, n = ncol(m))
}

log_msg("Reading ", nrow(samples), " GSE216783 samples (raw 10x, cell-calling)")
mats <- list(); meta_rows <- list()
for (i in seq_len(nrow(samples))) {
  r <- read_qc_sample(samples$gsm[i], samples$sample[i], samples$assay_type[i])
  if (is.null(r) || r$n < 50) next
  mats[[samples$sample[i]]] <- r$mat
  meta_rows[[samples$sample[i]]] <- data.frame(
    barcode = colnames(r$mat), sample = samples$sample[i],
    assay_type = samples$assay_type[i], stringsAsFactors = FALSE)
  invisible(gc())
}
if (length(mats) < 2) stop("Too few samples passed QC")

# 共同基因合并
common_genes <- Reduce(intersect, lapply(mats, rownames))
log_msg("Common genes across samples: ", length(common_genes))
mats <- lapply(mats, function(x) x[common_genes, , drop = FALSE])
counts <- do.call(cbind, mats)
meta <- do.call(rbind, meta_rows)
rownames(meta) <- meta$barcode
rm(mats, meta_rows); invisible(gc())
log_msg("Merged matrix: ", nrow(counts), " genes x ", ncol(counts), " cells")

# ---- 归一化 (log1p CPM/1e4) ----
libsize <- Matrix::colSums(counts)
norm <- counts
norm@x <- log1p(norm@x / rep(libsize, diff(norm@p)) * 1e4)

# ---- marker-score per cell → 注释 (argmax; 需 top score>0) ----
score_mat <- sapply(CELLTYPES, function(ct) {
  g <- intersect(MARKERS[[ct]], rownames(norm))
  if (!length(g)) return(rep(0, ncol(norm)))
  Matrix::colMeans(norm[g, , drop = FALSE])
})
rownames(score_mat) <- colnames(norm)
top <- max.col(score_mat, ties.method = "first")
topval <- score_mat[cbind(seq_len(nrow(score_mat)), top)]
meta$cell_type <- ifelse(topval > 0, CELLTYPES[top], "Unassigned")
log_msg("Compartment assignment done")
print(table(meta$cell_type, meta$assay_type))

write.csv(meta, file.path(outdir, "cell_meta.csv"), row.names = FALSE)
comp <- as.data.frame(table(cell_type = meta$cell_type))
comp$percent <- round(100 * comp$Freq / sum(comp$Freq), 2)
write.csv(comp, file.path(outdir, "compartment_composition.csv"), row.names = FALSE)

# ---- 区室平均表达参考 (供反卷积) ----
real_cells <- meta$cell_type %in% CELLTYPES
meanexpr <- sapply(CELLTYPES, function(ct) {
  cols <- which(meta$cell_type == ct)
  if (!length(cols)) return(rep(NA_real_, nrow(norm)))
  Matrix::rowMeans(norm[, cols, drop = FALSE])
})
rownames(meanexpr) <- rownames(norm)
write.csv(data.frame(gene = rownames(meanexpr), meanexpr, check.names = FALSE),
          file.path(outdir, "celltype_meanexpr_gse216783.csv"), row.names = FALSE)

# ---- 候选基因区室定位 + 与 GSE230375 比较 ----
cand_present <- intersect(CANDIDATES, rownames(norm))
loc <- do.call(rbind, lapply(cand_present, function(g) {
  vals <- sapply(CELLTYPES, function(ct) {
    cols <- which(meta$cell_type == ct)
    if (!length(cols)) return(c(mean = NA, pct = NA))
    c(mean = mean(norm[g, cols]), pct = 100 * mean(counts[g, cols] > 0))
  })
  dom <- CELLTYPES[which.max(vals["mean", ])]
  data.frame(gene = g, dominant_celltype_gse216783 = dom,
             dom_mean = round(max(vals["mean", ], na.rm = TRUE), 3),
             stringsAsFactors = FALSE)
}))
# GSE230375 主导细胞 (来自加深1的 meanexpr)
ref230 <- read.csv(file.path(procd, "GSE230375_celltype_meanexpr.csv"), check.names = FALSE)
ref230_map <- c("Schwann/肿瘤细胞"="Schwann/tumour","巨噬/髓系"="Myeloid","T细胞"="T cell",
                "成纤维"="Fibroblast","内皮"="Endothelial","周细胞"="Pericyte")
r230 <- ref230[match(loc$gene, ref230$gene), ]
ct_cols <- intersect(names(ref230_map), colnames(ref230))
dom230 <- apply(r230[, ct_cols, drop = FALSE], 1, function(z) {
  if (all(is.na(z))) return(NA_character_)
  ref230_map[ct_cols[which.max(z)]]
})
loc$dominant_celltype_gse230375 <- unname(dom230)
loc$concordant <- loc$dominant_celltype_gse216783 == loc$dominant_celltype_gse230375
write.csv(loc, file.path(outdir, "candidate_compartment_localization.csv"), row.names = FALSE)
log_msg("Candidate localization concordance: ",
        sum(loc$concordant, na.rm = TRUE), "/", nrow(loc))

# ---- 第二参考 NNLS 反卷积 (敏感性) ----
to_linear <- function(x) { x <- as.matrix(x); 2^(x - min(x, na.rm = TRUE)) }
select_markers <- function(mn, n = 200) {
  lin <- expm1(pmax(as.matrix(mn), 0))
  mx <- sapply(seq_len(ncol(lin)), function(j) apply(lin[, -j, drop = FALSE], 1, max))
  spec <- log2((lin + 0.05) / (mx + 0.05))
  unique(unlist(lapply(seq_len(ncol(lin)), function(j) {
    el <- which(lin[, j] > 0.05 & is.finite(spec[, j]))
    head(rownames(lin)[el[order(spec[el, j], decreasing = TRUE)]], n)
  })))
}
run_nnls <- function(ref, bulk) {
  common <- intersect(rownames(ref), rownames(bulk))
  ref <- as.matrix(ref[common, ]); bul <- as.matrix(bulk[common, ])
  keep <- rowSums(ref) > 0 & apply(bul, 1, function(z) all(is.finite(z)))
  ref <- ref[keep, ]; bul <- bul[keep, ]
  est <- t(vapply(seq_len(ncol(bul)), function(i) {
    p <- pmax(stats::coef(nnls::nnls(ref, bul[, i])), 0)
    if (!sum(p) > 0) p <- rep(1, length(p)); p / sum(p)
  }, numeric(ncol(ref))))
  rownames(est) <- colnames(bul); colnames(est) <- colnames(ref); est
}
read_expr <- function(p) { x <- read.csv(p, check.names = FALSE); m <- as.matrix(x[,-1]); rownames(m) <- x[[1]]; m }

suppressMessages(library(nnls))
mk <- select_markers(meanexpr)
ref216 <- expm1(pmax(meanexpr[intersect(mk, rownames(meanexpr)), , drop = FALSE], 0))

cohorts <- list(GSE39645 = "GSE39645_expr_gene.csv", GSE108524 = "GSE108524_expr_gene.csv")
cmp_rows <- list()
for (co in names(cohorts)) {
  expr <- read_expr(file.path(procd, cohorts[[co]]))
  lin <- to_linear(expr)
  est <- run_nnls(ref216, lin)
  write.csv(data.frame(sample = rownames(est), est, check.names = FALSE),
            file.path(outdir, paste0("cell_fractions_", co, "_gse216783ref.csv")), row.names = FALSE)
  # 与 GSE230375 参考 NNLS (13 产物) 比较 sample-level Spearman per celltype
  old_f <- file.path(root, "results", "deconv", paste0("cell_fractions_", co, "_nnls.csv"))
  if (file.exists(old_f)) {
    old <- read.csv(old_f, check.names = FALSE)
    map230 <- c("Schwann/肿瘤细胞"="Schwann/tumour","巨噬/髓系"="Myeloid","T细胞"="T cell",
                "成纤维"="Fibroblast","内皮"="Endothelial","周细胞"="Pericyte")
    for (ctn in colnames(est)) {
      oldcol <- names(map230)[map230 == ctn]
      if (length(oldcol) && oldcol %in% colnames(old)) {
        m <- merge(data.frame(sample = rownames(est), new = est[, ctn]),
                   old[, c("sample", oldcol)], by = "sample")
        r <- suppressWarnings(cor(m$new, m[[oldcol]], method = "spearman"))
        cmp_rows[[paste(co, ctn)]] <- data.frame(cohort = co, cell_type = ctn,
          n = nrow(m), spearman_vs_GSE230375ref = round(r, 3))
      }
    }
  }
}
cmp <- do.call(rbind, cmp_rows)
write.csv(cmp, file.path(outdir, "reference_comparison_deconv.csv"), row.names = FALSE)

# ---- 报告 ----
fmt <- function(x) {
  if (is.null(x) || !nrow(x)) return("_无_")
  paste(c(paste0("| ", paste(names(x), collapse=" | "), " |"),
          paste0("|", paste(rep("---", ncol(x)), collapse="|"), "|"),
          apply(x,1,function(z) paste0("| ", paste(z, collapse=" | "), " |"))), collapse="\n")
}
report <- c(
  "# GSE216783 second-reference / within-VS localization report (TASK-D)",
  "",
  "## Scope & wording",
  "GSE216783 (Barrett et al. 2024, doi:10.1038/s41467-023-42762-w; 15 sporadic VS, scRNA+snRNA, Cell Ranger RAW matrices) 仅用于 (1) VS 内部独立的候选区室定位; (2) 第二个 tumour-derived 参考的反卷积敏感性分析。**无正常神经对照, 且与已引 Barrett 同源**, 故不作独立 VS-vs-nerve 比较、不作独立生物学验证、不宣称破除单参考循环。证据为 in-silico, 措辞 候选/与…一致/提示。",
  "",
  paste0("## QC & 区室构成 (cell-calling: nCount>=", MIN_COUNTS, ", nFeature>=", MIN_FEAT,
         ", %MT<", MAX_MT_PCT, ")"),
  paste0("- 合并后细胞数: ", nrow(meta), " (", length(unique(meta$sample)), " 样本; ",
         sum(meta$assay_type=="scRNA"), " scRNA + ", sum(meta$assay_type=="snRNA"), " snRNA cells)。"),
  "",
  fmt(comp),
  "",
  "## 候选基因区室定位: GSE216783 vs GSE230375",
  paste0("一致 ", sum(loc$concordant, na.rm=TRUE), "/", nrow(loc), " 个候选的主导区室在两套独立 scRNA 间相同。"),
  "",
  fmt(loc),
  "",
  "## 第二参考反卷积敏感性 (GSE216783-ref NNLS vs GSE230375-ref NNLS, sample-level Spearman)",
  fmt(cmp),
  "",
  "## 产物 (绝对路径)",
  paste0("- 脚本: `", file.path(root,"scripts","14_gse216783_second_ref.R"), "`"),
  paste0("- 细胞 meta: `", file.path(outdir,"cell_meta.csv"), "`"),
  paste0("- 区室构成: `", file.path(outdir,"compartment_composition.csv"), "`"),
  paste0("- 候选定位: `", file.path(outdir,"candidate_compartment_localization.csv"), "`"),
  paste0("- 反卷积敏感性: `", file.path(outdir,"reference_comparison_deconv.csv"), "`"),
  paste0("- 参考均表达: `", file.path(outdir,"celltype_meanexpr_gse216783.csv"), "`"),
  paste0("- 原始数据: `", rawdir, "` (15 GSM, Cell Ranger RAW)"),
  "",
  "## 局限",
  "- GSE216783 全为 sporadic VS, 无正常神经/对照 → 不能独立复现 VS-vs-nerve 组成差异。",
  "- 与 Barrett 2024 同源 → 第二参考缓解算法依赖, 不构成生物学独立验证。",
  "- RAW 矩阵自做 cell-calling(阈值法, 非 emptyDrops); 粗注释用 marker-score argmax, 未做双细胞/ambient 校正; scRNA+snRNA 混合。",
  "- 反卷积为跨平台(scRNA 参考→microarray)探索性应用, 与加深2 同样限制。"
)
writeLines(report, file.path(docdir, "gse216783_report.md"), useBytes = TRUE)
log_msg("TASK-D 完成; 报告: ", file.path(docdir, "gse216783_report.md"))
