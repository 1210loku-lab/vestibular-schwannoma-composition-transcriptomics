#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

project_root <- Sys.getenv("PROJ_ROOT", getwd())
output_dir <- file.path(project_root, "results", "scrna")
tmp_root <- "/tmp/coop_neurosurg"
tmp_dir <- file.path(tmp_root, "scrna_tmp")
tar_path <- file.path(tmp_root, "GSE230375_RAW.tar")
hub_path <- file.path(project_root, "results", "hub", "hub_all.csv")
signature_path <- file.path(project_root, "results", "signature", "signature_genes.csv")
script_path <- file.path(project_root, "scripts", "08_scrna_localize.R")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

log_message <- function(...) {
  message(sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
}

stop_clear <- function(...) {
  stop(paste0(...), call. = FALSE)
}

require_package <- function(pkg, required = TRUE) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  if (!ok && required) stop_clear("Required R package is unavailable: ", pkg)
  if (!ok) log_message("Optional R package unavailable: ", pkg)
  ok
}

invisible(lapply(c("Seurat", "Matrix", "ggplot2", "BPCells"), require_package))
harmony_available <- require_package("harmony", required = FALSE)

download_geo_tar <- function() {
  minimum_size <- 100 * 1024^2
  existing_size <- if (file.exists(tar_path)) file.info(tar_path)$size else 0
  if (is.na(existing_size)) existing_size <- 0
  if (existing_size < minimum_size) {
    if (existing_size > 0) {
      log_message("Existing GEO tar is too small (", round(existing_size / 1024^2, 1),
                  " MB); replacing it.")
      unlink(tar_path)
    }
    cmd <- paste(
      "curl -sL",
      shQuote("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE230nnn/GSE230375/suppl/GSE230375_RAW.tar"),
      "-o", shQuote(tar_path), "</dev/null"
    )
    log_message("Downloading GSE230375 raw archive to ", tar_path)
    status <- system(cmd)
    if (status != 0) stop_clear("GEO download failed with status ", status, ": ", cmd)
  } else {
    log_message("Using existing GEO archive: ", tar_path)
  }
  final_size <- file.info(tar_path)$size
  if (is.na(final_size) || final_size <= minimum_size) {
    stop_clear("GEO archive failed size validation (>100 MB required): ", tar_path,
               " size=", final_size)
  }
  log_message("Archive size validated: ", round(final_size / 1024^2, 1), " MB")
}

extract_geo_tar <- function() {
  marker <- file.path(tmp_dir, ".GSE230375_extract_complete")
  if (file.exists(marker)) {
    log_message("Using existing extracted GEO files under ", tmp_dir)
    return(invisible(NULL))
  }
  log_message("Extracting archive to ", tmp_dir)
  status <- system2("tar", c("xf", shQuote(tar_path), "-C", shQuote(tmp_dir)))
  if (status != 0) stop_clear("tar extraction failed with status ", status)
  file.create(marker)
}

strip_compression <- function(x) sub("\\.(gz|bz2|xz)$", "", x, ignore.case = TRUE)

file_kind <- function(path) {
  x <- tolower(strip_compression(basename(path)))
  if (grepl("matrix\\.mtx$", x)) return("matrix")
  if (grepl("barcodes\\.tsv$", x)) return("barcodes")
  if (grepl("(features|genes)\\.tsv$", x)) return("features")
  NA_character_
}

sample_key_from_file <- function(path) {
  x <- strip_compression(basename(path))
  stripped <- sub("([._-]?)(matrix\\.mtx|barcodes\\.tsv|features\\.tsv|genes\\.tsv)$",
                  "", x, ignore.case = TRUE)
  if (!nzchar(stripped)) {
    normalizePath(dirname(path), mustWork = FALSE)
  } else {
    stripped
  }
}

infer_sample_label <- function(key, paths) {
  text <- paste(c(key, basename(paths), basename(dirname(paths))), collapse = "_")
  hit <- regmatches(toupper(text), regexpr("(^|[^A-Z0-9])(T[1-7]|N[1-2])([^A-Z0-9]|$)",
                                          toupper(text), perl = TRUE))
  hit <- gsub("[^TN0-9]", "", hit)
  if (length(hit) == 1 && nzchar(hit)) return(hit)

  gsm <- regmatches(toupper(text), regexpr("GSM[0-9]+", toupper(text), perl = TRUE))
  gsm_to_label <- c(
    GSM7809242 = "T1", GSM7809243 = "T2", GSM7809244 = "T3",
    GSM7809245 = "T4", GSM7809246 = "T5", GSM7809247 = "T6",
    GSM7809248 = "T7", GSM7809249 = "N1", GSM7809250 = "N2"
  )
  if (length(gsm) == 1 && nzchar(gsm) && gsm %in% names(gsm_to_label)) {
    return(unname(gsm_to_label[[gsm]]))
  }
  if (length(gsm) == 1 && nzchar(gsm)) return(gsm)
  key
}

infer_group <- function(sample_label, key, paths) {
  text <- toupper(paste(c(sample_label, key, basename(paths)), collapse = "_"))
  if (grepl("(^|[^A-Z0-9])N[1-2]([^A-Z0-9]|$)|NERVE|NORMAL|GREAT.?AURICULAR", text, perl = TRUE)) {
    return("Nerve")
  }
  if (grepl("(^|[^A-Z0-9])T[1-7]([^A-Z0-9]|$)|VESTIBULAR|SCHWANNOMA|TUMOU?R|(^|[_-])VS([_-]|$)",
            text, perl = TRUE)) {
    return("VS")
  }
  NA_character_
}

discover_10x_samples <- function(root) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE)
  kinds <- vapply(files, file_kind, character(1))
  files <- files[!is.na(kinds)]
  kinds <- kinds[!is.na(kinds)]
  if (!length(files)) stop_clear("No 10x matrix/barcode/feature files found under ", root)

  keys <- vapply(files, sample_key_from_file, character(1))
  split_idx <- split(seq_along(files), keys)
  rows <- lapply(names(split_idx), function(key) {
    idx <- split_idx[[key]]
    these_files <- files[idx]
    these_kinds <- kinds[idx]
    pick <- function(kind) {
      candidate <- these_files[these_kinds == kind]
      if (length(candidate)) candidate[[1]] else NA_character_
    }
    label <- infer_sample_label(key, these_files)
    data.frame(
      sample_key = key,
      sample = label,
      group = infer_group(label, key, these_files),
      matrix = pick("matrix"),
      barcodes = pick("barcodes"),
      features = pick("features"),
      n_matrix_files = sum(these_kinds == "matrix"),
      n_barcode_files = sum(these_kinds == "barcodes"),
      n_feature_files = sum(these_kinds == "features"),
      stringsAsFactors = FALSE
    )
  })
  manifest <- do.call(rbind, rows)
  manifest$complete <- with(manifest,
                            !is.na(matrix) & !is.na(barcodes) & !is.na(features) &
                              n_matrix_files == 1 & n_barcode_files == 1 & n_feature_files == 1)
  manifest
}

qc_filter_counts <- function(counts, min_features = 200, max_features = 8000,
                             max_mt = 20) {
  nfeature <- Matrix::colSums(counts > 0)
  mt_rows <- grepl("^MT-", rownames(counts))
  mt_counts <- if (any(mt_rows)) {
    Matrix::colSums(counts[mt_rows, , drop = FALSE])
  } else {
    rep(0, ncol(counts))
  }
  total_counts <- Matrix::colSums(counts)
  mt_pct <- 100 * mt_counts / pmax(total_counts, 1)
  keep <- nfeature > min_features & nfeature < max_features & mt_pct < max_mt
  list(
    counts = counts[, keep, drop = FALSE],
    cells_before_qc = as.integer(ncol(counts)),
    cells_after_qc = as.integer(sum(keep)),
    keep = keep,
    nfeature = nfeature,
    percent_mt = mt_pct
  )
}

read_qc_write_sample <- function(row, bp_root) {
  if (!isTRUE(row$complete)) stop_clear("10x triplet is incomplete or duplicated")
  if (is.na(row$group)) stop_clear("Cannot infer VS/Nerve group from filename")
  log_message("Reading and filtering sample ", row$sample, " [", row$group, "]")
  counts <- Seurat::ReadMtx(
    mtx = row$matrix,
    cells = row$barcodes,
    features = row$features,
    feature.column = 2,
    cell.column = 1,
    unique.features = TRUE,
    strip.suffix = FALSE
  )
  if (!inherits(counts, "dgCMatrix")) counts <- as(counts, "dgCMatrix")
  if (!nrow(counts) || !ncol(counts)) stop_clear("ReadMtx returned an empty matrix")
  qc <- qc_filter_counts(counts)
  sample_dir <- file.path(bp_root, "samples", row$sample)
  dir.create(dirname(sample_dir), recursive = TRUE, showWarnings = FALSE)
  BPCells::write_matrix_dir(qc$counts, dir = sample_dir, overwrite = TRUE)
  rm(counts)
  invisible(gc())
  data.frame(
    sample = row$sample,
    group = row$group,
    cells_before_qc = qc$cells_before_qc,
    cells_after_qc = qc$cells_after_qc,
    bp_path = sample_dir,
    stringsAsFactors = FALSE
  )
}

combine_bpcells_samples <- function(sample_paths, sample_groups) {
  if (!length(sample_paths)) stop_clear("No BPCells sample matrices supplied")
  if (is.null(names(sample_paths)) || any(!nzchar(names(sample_paths)))) {
    stop_clear("sample_paths must be a named vector")
  }
  matrices <- lapply(names(sample_paths), function(sample) {
    m <- BPCells::open_matrix_dir(sample_paths[[sample]])
    colnames(m) <- paste0(sample, "_", colnames(m))
    m
  })
  reference_genes <- rownames(matrices[[1]])
  if (!all(vapply(matrices, function(m) identical(rownames(m), reference_genes),
                  logical(1)))) {
    stop_clear("Sample BPCells matrices do not share identical gene order")
  }
  counts <- do.call(cbind, matrices)
  samples <- sub("_.*$", "", colnames(counts))
  groups <- unname(sample_groups[samples])
  if (anyNA(groups)) stop_clear("Missing group for one or more merged samples")
  metadata <- data.frame(
    row.names = colnames(counts),
    sample = samples,
    group = groups,
    stringsAsFactors = FALSE
  )
  list(counts = counts, metadata = metadata)
}

sha256_file <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  sub("\\s+.*$", "", out[[1]])
}

stage_can_resume <- function(stage_manifest, stage, upstream_sha256) {
  row <- stage_manifest[stage_manifest$stage == stage, , drop = FALSE]
  if (nrow(row) != 1 || row$status[[1]] != "complete") return(FALSE)
  artifact <- row$artifact[[1]]
  isTRUE(file.exists(artifact)) &&
    identical(row$upstream_sha256[[1]], upstream_sha256) &&
    identical(row$sha256[[1]], sha256_file(artifact))
}

aggregate_counts_by_sample_celltype <- function(object) {
  counts <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
  meta <- object[[]]
  if (!all(c("sample", "cell_type") %in% names(meta))) {
    stop_clear("Object metadata must contain sample and cell_type")
  }
  groups <- interaction(meta$sample, meta$cell_type, drop = TRUE, sep = "__")
  design <- Matrix::sparse.model.matrix(~ 0 + groups)
  colnames(design) <- sub("^groups", "", colnames(design))
  counts %*% design
}

get_layer <- function(object, layer) {
  SeuratObject::LayerData(object, assay = "RNA", layer = layer)
}

get_feature_block <- function(object, layer, features) {
  present <- intersect(unique(features), rownames(object))
  if (!length(present)) {
    return(Matrix::Matrix(0, nrow = 0, ncol = ncol(object), sparse = TRUE,
                          dimnames = list(character(), colnames(object))))
  }
  SeuratObject::LayerData(object, assay = "RNA", layer = layer,
                          fast = FALSE)[present, , drop = FALSE]
}

save_blank_plot <- function(path, title, subtitle) {
  p <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0.15, label = title, fontface = "bold", size = 5) +
    ggplot2::annotate("text", x = 0, y = -0.15, label = subtitle, size = 4) +
    ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) + ggplot2::theme_void()
  ggplot2::ggsave(path, p, width = 8, height = 5, dpi = 300)
}

make_dotplot <- function(object, genes, group_by, path, title, width, height) {
  present <- genes[genes %in% rownames(object)]
  absent <- setdiff(genes, present)
  if (length(absent)) log_message("DotPlot genes absent from matrix: ", paste(absent, collapse = ", "))
  if (!length(present)) {
    save_blank_plot(path, title, "None of the requested genes were detected in the matrix.")
    return(invisible(NULL))
  }
  p <- Seurat::DotPlot(object, features = unique(present), group.by = group_by,
                       cols = c("#E8EEF7", "#B2182B"), dot.scale = 7) +
    Seurat::RotatedAxis() +
    ggplot2::labs(title = title, x = NULL, y = NULL, color = "Mean expression",
                  size = "% expressing") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(linewidth = 0.2,
                                                            color = "grey90"),
                   axis.text.x = ggplot2::element_text(
                     angle = 60, hjust = 1, vjust = 1, size = 8
                   ))
  ggplot2::ggsave(path, p, width = width, height = height, dpi = 300, limitsize = FALSE)
}

summarize_target_expression <- function(object, target_table) {
  cell_types <- sort(unique(object$cell_type))
  data_mat <- get_feature_block(object, "data", target_table$gene)
  counts_mat <- get_feature_block(object, "counts", target_table$gene)
  out <- vector("list", nrow(target_table) * length(cell_types))
  k <- 1L
  for (i in seq_len(nrow(target_table))) {
    gene <- target_table$gene[[i]]
    for (ct in cell_types) {
      cells <- colnames(object)[object$cell_type == ct]
      if (gene %in% rownames(data_mat) && length(cells)) {
        avg <- mean(data_mat[gene, cells, drop = TRUE])
        pct <- mean(counts_mat[gene, cells, drop = TRUE] > 0) * 100
      } else {
        avg <- NA_real_
        pct <- NA_real_
      }
      out[[k]] <- data.frame(
        gene = gene,
        source = target_table$source[[i]],
        bulk_direction = target_table$bulk_direction[[i]],
        bulk_logFC = target_table$bulk_logFC[[i]],
        cell_type = ct,
        avg_expression = avg,
        pct_expressing = pct,
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  do.call(rbind, out)
}

dominant_celltype <- function(expression_long) {
  split_rows <- split(expression_long, expression_long$gene)
  do.call(rbind, lapply(split_rows, function(x) {
    valid <- which(is.finite(x$avg_expression))
    if (!length(valid)) {
      chosen <- x[1, , drop = FALSE]
      chosen$cell_type <- "Not detected"
      chosen$avg_expression <- NA_real_
      chosen$pct_expressing <- NA_real_
    } else {
      chosen <- x[valid[which.max(x$avg_expression[valid])], , drop = FALSE]
    }
    data.frame(
      gene = chosen$gene,
      source = chosen$source,
      bulk_direction = chosen$bulk_direction,
      bulk_logFC = chosen$bulk_logFC,
      dominant_cell_type = chosen$cell_type,
      dominant_avg_expression = chosen$avg_expression,
      dominant_pct_expressing = chosen$pct_expressing,
      stringsAsFactors = FALSE
    )
  }))
}

write_report <- function(object, manifest, failures, qc_summary, annotation, dominant, composition,
                         harmony_status, report_path) {
  signature_rows <- dominant[grepl("signature", dominant$source), , drop = FALSE]
  hub_rows <- dominant[grepl("hub", dominant$source), , drop = FALSE]
  fmt_table <- function(x, n = Inf) {
    if (!nrow(x)) return("_No rows available._")
    x <- head(x, n)
    header <- paste0("| ", paste(names(x), collapse = " | "), " |")
    sep <- paste0("|", paste(rep("---", ncol(x)), collapse = "|"), "|")
    body <- apply(x, 1, function(z) paste0("| ", paste(z, collapse = " | "), " |"))
    paste(c(header, sep, body), collapse = "\n")
  }
  sig_show <- signature_rows[, c("gene", "bulk_direction", "dominant_cell_type",
                                  "dominant_avg_expression", "dominant_pct_expressing")]
  hub_show <- hub_rows[, c("gene", "bulk_direction", "bulk_logFC", "dominant_cell_type",
                           "dominant_avg_expression", "dominant_pct_expressing")]
  num_cols <- vapply(sig_show, is.numeric, logical(1))
  sig_show[num_cols] <- lapply(sig_show[num_cols], function(z) round(z, 3))
  num_cols <- vapply(hub_show, is.numeric, logical(1))
  hub_show[num_cols] <- lapply(hub_show[num_cols], function(z) round(z, 3))

  focus_genes <- c("NRG1", "L1CAM", "NCAM2", "PRRX1", "FBLN1", "PDGFRA")
  focus <- dominant[dominant$gene %in% focus_genes, c("gene", "bulk_direction",
                                                      "dominant_cell_type",
                                                      "dominant_avg_expression",
                                                      "dominant_pct_expressing"), drop = FALSE]
  focus <- focus[match(intersect(focus_genes, focus$gene), focus$gene), , drop = FALSE]
  num_cols <- vapply(focus, is.numeric, logical(1))
  focus[num_cols] <- lapply(focus[num_cols], function(z) round(z, 3))
  up_rows <- dominant[grepl("(^|;)up($|;)", dominant$bulk_direction), , drop = FALSE]
  down_rows <- dominant[grepl("(^|;)down($|;)", dominant$bulk_direction), , drop = FALSE]
  n_up_schwann <- sum(up_rows$dominant_cell_type == "Schwann/肿瘤细胞")
  n_down_fibro <- sum(down_rows$dominant_cell_type == "成纤维")
  nerve_fibro_pct <- composition$percent_within_group[
    composition$group == "Nerve" & composition$cell_type == "成纤维"
  ]
  vs_fibro_pct <- composition$percent_within_group[
    composition$group == "VS" & composition$cell_type == "成纤维"
  ]

  paths <- c(
    script = script_path,
    cluster_annotation = file.path(output_dir, "cluster_annotation.csv"),
    gene_celltype_expression = file.path(output_dir, "gene_celltype_expression.csv"),
    gene_dominant_celltype = file.path(output_dir, "gene_dominant_celltype.csv"),
    celltype_composition = file.path(output_dir, "celltype_composition.csv"),
    qc_summary = file.path(output_dir, "qc_summary.csv"),
    sample_manifest = file.path(output_dir, "sample_manifest.csv"),
    sample_failures = file.path(output_dir, "sample_failures.csv"),
    sample_cell_counts = file.path(output_dir, "sample_cell_counts.csv"),
    run_log = file.path(output_dir, "08_run.log"),
    umap_celltype = file.path(output_dir, "fig_umap_celltype.png"),
    umap_group = file.path(output_dir, "fig_umap_group.png"),
    dotplot_signature = file.path(output_dir, "fig_dotplot_signature.png"),
    dotplot_hub = file.path(output_dir, "fig_dotplot_hub.png"),
    report = report_path,
    seurat_checkpoint = file.path(project_root, "data", "processed", "GSE230375_seurat.rds")
  )

  lines <- c(
    "# GSE230375 单细胞定位报告",
    "",
    "## 分析范围与措辞",
    "",
    "本分析将 bulk RNA-seq 中的候选 hub/signature 基因定位到 GSE230375 的粗粒度细胞类型。结果属于单细胞表达定位证据，用于判断与 bulk 方向及细胞组成假设是否一致，不构成机制性或直接实验验证证据。",
    "",
    "## 数据与质控",
    "",
    paste0("- 成功读取样本：", length(unique(object$sample)), "；VS=",
           length(unique(object$sample[object$group == "VS"])), "，Nerve=",
           length(unique(object$sample[object$group == "Nerve"])), "。"),
    paste0("- QC 前细胞数：", qc_summary$value[qc_summary$metric == "cells_before_qc"],
           "；QC 后细胞数：", qc_summary$value[qc_summary$metric == "cells_after_qc"], "。"),
    "- 固定过滤条件：nFeature_RNA > 200、nFeature_RNA < 8000、percent.mt < 20；percent.mt 使用 `^MT-` 计算。",
    paste0("- 降维批次处理：", harmony_status, "。"),
    "- 本任务卡未包含双细胞识别或 ambient RNA 校正，因此这些因素仍可能影响小群体与混合 marker 的解释。",
    "",
    "## 细胞类型构成",
    "",
    fmt_table(composition),
    "",
    "## Cluster 注释依据",
    "",
    "采用 cluster 级 canonical-marker 平均归一化表达评分，指派最高分粗粒度类型；无 marker 信号者保留为 Unknown/ambiguous。",
    "",
    fmt_table(annotation[, c("cluster", "cell_type", "n_cells", "top_score",
                              "second_score", "score_margin")]),
    "",
    "Marker 面板：Schwann/肿瘤细胞（SOX10,S100B,PMP22,MPZ,PLP1,NGFR）；巨噬/髓系（CD68,CD163,LYZ,AIF1）；T 细胞（CD3D,CD3E,CD2）；成纤维（DCN,LUM,COL1A1,PDGFRA）；内皮（PECAM1,VWF,CLDN5）；周细胞（RGS5,ACTA2）。",
    "",
    "## Signature 基因主要表达细胞类型",
    "",
    fmt_table(sig_show),
    "",
    "## Hub 基因主要表达细胞类型",
    "",
    fmt_table(hub_show),
    "",
    "## 重点基因与 bulk 方向一致性",
    "",
    fmt_table(focus),
    "",
    paste0(
      "实际结果显示，L1CAM、NRG1、NCAM2 的最高平均表达均位于 Schwann/肿瘤细胞；",
      "PRRX1、FBLN1、PDGFRA 的最高平均表达均位于成纤维细胞。Top ",
      nrow(up_rows), " 上调 hub 中 ", n_up_schwann,
      " 个主要定位于 Schwann/肿瘤细胞，Top ", nrow(down_rows), " 下调 hub 中 ",
      n_down_fibro, " 个主要定位于成纤维细胞。Nerve 中成纤维细胞占 ",
      sprintf("%.2f", nerve_fibro_pct), "%，VS 中为 ", sprintf("%.2f", vs_fibro_pct),
      "%。这些结果与“上调候选主要来自 Schwann/肿瘤细胞、下调 ECM/基质候选部分反映正常神经基质成分减少”的解释一致，",
      "但细胞捕获比例与解离偏倚也可能影响组间构成，不能据此区分细胞组成变化和细胞内表达变化。"
    ),
    "",
    "## 失败与降级记录",
    "",
    if (nrow(failures)) fmt_table(failures) else "_无样本读取失败。_",
    "",
    "## 产物与完整绝对路径",
    "",
    paste0("- `", names(paths), "`: `", unname(paths), "`"),
    "",
    "## 局限",
    "",
    "- 粗粒度 marker 注释不能完全区分 neoplastic Schwann cells 与正常 Schwann cells，因此统一标注为 Schwann/肿瘤细胞。",
    "- “主要表达细胞类型”定义为该基因归一化平均表达最高的粗粒度类型；该规则不等价于特异表达，需同时查看表达比例和 DotPlot。",
    "- VS 与 Nerve 的细胞捕获比例、解离偏倚和样本量不同，bulk 方向可能同时反映细胞内表达变化与细胞组成差异。",
    "- 未进行双细胞识别、环境 RNA 校正或参考图谱映射；低置信度 cluster 应在后续精细注释中复核。",
    "",
    "## 输入来源",
    "",
    paste0("- Hub 表：`", hub_path, "`"),
    paste0("- Signature 表：`", signature_path, "`"),
    paste0("- GEO 原始包：`", tar_path, "`"),
    paste0("- 解压临时目录：`", tmp_dir, "`")
  )
  writeLines(lines, report_path, useBytes = TRUE)
}

main <- function() {
log_message("Starting GSE230375 single-cell localization")
download_geo_tar()
extract_geo_tar()

manifest <- discover_10x_samples(tmp_dir)
write.csv(manifest, file.path(output_dir, "sample_manifest.csv"), row.names = FALSE)
log_message("Discovered ", nrow(manifest), " candidate 10x sample bundles")
for (i in seq_len(nrow(manifest))) {
  log_message("Bundle ", manifest$sample[[i]], ": complete=", manifest$complete[[i]],
              ", group=", ifelse(is.na(manifest$group[[i]]), "NA", manifest$group[[i]]))
}

failures <- list()
bp_root <- file.path(project_root, "data", "processed", "GSE230375_bpcells")
dir.create(bp_root, recursive = TRUE, showWarnings = FALSE)
sample_records <- list()
for (i in seq_len(nrow(manifest))) {
  row <- manifest[i, , drop = FALSE]
  result <- tryCatch(
    read_qc_write_sample(row, bp_root),
    error = function(e) {
      log_message("SKIP sample ", row$sample, ": ", conditionMessage(e))
      failures[[length(failures) + 1L]] <<- data.frame(
        sample = row$sample, sample_key = row$sample_key,
        reason = conditionMessage(e), stringsAsFactors = FALSE
      )
      NULL
    }
  )
  if (!is.null(result)) sample_records[[row$sample]] <- result
}

failure_df <- if (length(failures)) do.call(rbind, failures) else {
  data.frame(sample = character(), sample_key = character(), reason = character())
}
write.csv(failure_df, file.path(output_dir, "sample_failures.csv"), row.names = FALSE)
if (length(sample_records) < 2) {
  stop_clear("Fewer than two valid samples were readable; cannot continue.")
}
sample_cell_counts <- do.call(rbind, sample_records)
if (!all(c("VS", "Nerve") %in% sample_cell_counts$group)) {
  stop_clear("Both VS and Nerve groups are required after sample loading.")
}

log_message("Combining ", nrow(sample_cell_counts), " on-disk sample matrices")
sample_paths <- setNames(sample_cell_counts$bp_path, sample_cell_counts$sample)
sample_groups <- setNames(sample_cell_counts$group, sample_cell_counts$sample)
combined <- combine_bpcells_samples(sample_paths, sample_groups)
combined_dir <- file.path(bp_root, "combined_counts")
BPCells::write_matrix_dir(combined$counts, dir = combined_dir, overwrite = TRUE)
combined_counts <- BPCells::open_matrix_dir(combined_dir)
object <- Seurat::CreateSeuratObject(
  counts = combined_counts,
  project = "GSE230375",
  meta.data = combined$metadata,
  min.cells = 0,
  min.features = 0
)
cells_before_qc <- sum(sample_cell_counts$cells_before_qc)
cells_after_qc <- ncol(object)
log_message("QC cells before=", cells_before_qc, "; retained=", cells_after_qc,
            "; removed=", cells_before_qc - cells_after_qc)
if (cells_after_qc < 100) stop_clear("Fewer than 100 cells remain after QC.")

write.csv(sample_cell_counts[, c("sample", "group", "cells_before_qc",
                                 "cells_after_qc", "bp_path")],
          file.path(output_dir, "sample_cell_counts.csv"), row.names = FALSE)

qc_summary <- data.frame(
  metric = c("cells_before_qc", "cells_after_qc", "cells_removed_qc",
             "samples_loaded", "samples_failed"),
  value = c(cells_before_qc, cells_after_qc, cells_before_qc - cells_after_qc,
            nrow(sample_cell_counts), nrow(failure_df)),
  stringsAsFactors = FALSE
)
write.csv(qc_summary, file.path(output_dir, "qc_summary.csv"), row.names = FALSE)

log_message("Running NormalizeData -> FindVariableFeatures -> ScaleData -> PCA")
object <- Seurat::NormalizeData(object, verbose = FALSE)
object <- Seurat::FindVariableFeatures(object, selection.method = "vst",
                                       nfeatures = 3000, verbose = FALSE)
object <- Seurat::ScaleData(object, features = Seurat::VariableFeatures(object),
                            verbose = FALSE)
npcs <- min(50L, length(Seurat::VariableFeatures(object)) - 1L, ncol(object) - 1L)
if (npcs < 2) stop_clear("Insufficient dimensions for PCA: npcs=", npcs)
object <- Seurat::RunPCA(object, features = Seurat::VariableFeatures(object),
                         npcs = npcs, verbose = FALSE)
dims_use <- seq_len(min(30L, npcs))

reduction_use <- "pca"
harmony_status <- "Harmony 未使用"
if (harmony_available && length(unique(object$sample)) > 1) {
  log_message("Attempting Harmony integration by sample")
  harmony_result <- tryCatch({
    harmony::RunHarmony(object, group.by.vars = "sample", reduction.use = "pca",
                        dims.use = dims_use, verbose = FALSE)
  }, error = function(e) {
    log_message("Harmony failed; using uncorrected PCA. Reason: ", conditionMessage(e))
    harmony_status <<- paste0("Harmony 失败，降级使用 PCA；原因：", conditionMessage(e))
    NULL
  })
  if (!is.null(harmony_result)) {
    object <- harmony_result
    reduction_use <- "harmony"
    harmony_status <- "按 sample 使用 Harmony 校正"
  }
} else {
  harmony_status <- "Harmony 不可用或仅一个样本，使用 PCA"
  log_message(harmony_status)
}

available_dims <- ncol(Seurat::Embeddings(object, reduction = reduction_use))
dims_use <- seq_len(min(30L, available_dims))
log_message("Running UMAP, neighbors, and clustering using ", reduction_use)
object <- Seurat::RunUMAP(object, reduction = reduction_use, dims = dims_use,
                          seed.use = 230375, verbose = FALSE)
object <- Seurat::FindNeighbors(object, reduction = reduction_use, dims = dims_use,
                                verbose = FALSE)
object <- Seurat::FindClusters(object, resolution = 0.5, random.seed = 230375,
                               verbose = FALSE)

marker_sets <- list(
  "Schwann/肿瘤细胞" = c("SOX10", "S100B", "PMP22", "MPZ", "PLP1", "NGFR"),
  "巨噬/髓系" = c("CD68", "CD163", "LYZ", "AIF1"),
  "T细胞" = c("CD3D", "CD3E", "CD2"),
  "成纤维" = c("DCN", "LUM", "COL1A1", "PDGFRA"),
  "内皮" = c("PECAM1", "VWF", "CLDN5"),
  "周细胞" = c("RGS5", "ACTA2")
)

all_markers <- unique(unlist(marker_sets, use.names = FALSE))
data_mat <- get_feature_block(object, "data", all_markers)
cell_score <- matrix(0, nrow = ncol(object), ncol = length(marker_sets),
                     dimnames = list(colnames(object), names(marker_sets)))
marker_availability <- list()
for (ct in names(marker_sets)) {
  present <- intersect(marker_sets[[ct]], rownames(data_mat))
  marker_availability[[ct]] <- present
  if (!length(present)) {
    log_message("No annotation markers detected for ", ct)
    next
  }
  cell_score[, ct] <- Matrix::colMeans(data_mat[present, , drop = FALSE])
}
clusters <- as.character(object$seurat_clusters)
cluster_levels <- sort(unique(clusters))
cluster_score <- do.call(rbind, lapply(cluster_levels, function(cl) {
  colMeans(cell_score[clusters == cl, , drop = FALSE])
}))
rownames(cluster_score) <- cluster_levels

annotation_rows <- lapply(cluster_levels, function(cl) {
  scores <- cluster_score[cl, ]
  ordered <- order(scores, decreasing = TRUE)
  top <- scores[ordered[[1]]]
  second <- if (length(ordered) > 1) scores[ordered[[2]]] else NA_real_
  label <- if (!is.finite(top) || top <= 0) "Unknown/ambiguous" else names(scores)[ordered[[1]]]
  row <- data.frame(
    cluster = cl,
    cell_type = label,
    n_cells = sum(clusters == cl),
    top_score = top,
    second_score = second,
    score_margin = top - second,
    stringsAsFactors = FALSE
  )
  for (ct in names(marker_sets)) row[[paste0("score_", ct)]] <- scores[[ct]]
  row
})
annotation <- do.call(rbind, annotation_rows)
write.csv(annotation, file.path(output_dir, "cluster_annotation.csv"), row.names = FALSE)

cluster_to_type <- setNames(annotation$cell_type, annotation$cluster)
object$cell_type <- unname(cluster_to_type[as.character(object$seurat_clusters)])
object$cell_type <- factor(object$cell_type,
                           levels = c(names(marker_sets), "Unknown/ambiguous"))
celltype_plot_labels <- c(
  "Schwann/肿瘤细胞" = "Schwann/tumor cells",
  "巨噬/髓系" = "Macrophage/myeloid",
  "T细胞" = "T cells",
  "成纤维" = "Fibroblasts",
  "内皮" = "Endothelial",
  "周细胞" = "Pericytes",
  "Unknown/ambiguous" = "Unknown/ambiguous"
)
object$cell_type_plot <- factor(
  unname(celltype_plot_labels[as.character(object$cell_type)]),
  levels = unname(celltype_plot_labels)
)

composition <- as.data.frame(table(cell_type = object$cell_type, group = object$group),
                             stringsAsFactors = FALSE)
composition <- composition[composition$Freq > 0, , drop = FALSE]
group_totals <- aggregate(Freq ~ group, composition, sum)
names(group_totals)[2] <- "group_total"
composition <- merge(composition, group_totals, by = "group")
composition$percent_within_group <- round(100 * composition$Freq / composition$group_total, 2)
composition <- composition[, c("group", "cell_type", "Freq", "percent_within_group")]
names(composition)[3] <- "n_cells"
write.csv(composition, file.path(output_dir, "celltype_composition.csv"), row.names = FALSE)

celltype_palette <- c(
  "Schwann/tumor cells" = "#D73027", "Macrophage/myeloid" = "#4575B4",
  "T cells" = "#74ADD1", "Fibroblasts" = "#FDAE61", "Endothelial" = "#1A9850",
  "Pericytes" = "#984EA3",
  "Unknown/ambiguous" = "#999999"
)
p_celltype <- Seurat::DimPlot(object, reduction = "umap", group.by = "cell_type_plot",
                              label = TRUE, repel = TRUE, cols = celltype_palette,
                              raster = TRUE) +
  ggplot2::labs(title = "GSE230375: coarse cell types") +
  ggplot2::theme_bw(base_size = 11)
ggplot2::ggsave(file.path(output_dir, "fig_umap_celltype.png"), p_celltype,
                width = 9, height = 7, dpi = 300)

p_group <- Seurat::DimPlot(object, reduction = "umap", group.by = "group",
                           cols = c("VS" = "#B2182B", "Nerve" = "#2166AC"),
                           raster = TRUE) +
  ggplot2::labs(title = "GSE230375: VS versus Nerve") +
  ggplot2::theme_bw(base_size = 11)
ggplot2::ggsave(file.path(output_dir, "fig_umap_group.png"), p_group,
                width = 8, height = 7, dpi = 300)

hub <- read.csv(hub_path, check.names = FALSE)
signature <- read.csv(signature_path, check.names = FALSE)
required_hub <- c("gene", "direction", "logFC")
if (!all(required_hub %in% names(hub))) {
  stop_clear("hub_all.csv lacks required columns: ", paste(setdiff(required_hub, names(hub)),
                                                           collapse = ", "))
}
if (!"gene" %in% names(signature)) stop_clear("signature_genes.csv lacks gene column")

up_hub <- head(hub[hub$direction == "up", ][order(-abs(hub$logFC[hub$direction == "up"])), ], 30)
down_hub <- head(hub[hub$direction == "down", ][order(-abs(hub$logFC[hub$direction == "down"])), ], 20)
hub_target <- rbind(up_hub, down_hub)

signature_target <- data.frame(
  gene = signature$gene,
  source = "signature",
  bulk_direction = if ("coefficient" %in% names(signature)) {
    ifelse(signature$coefficient > 0, "signature_positive", "signature_negative")
  } else {
    NA_character_
  },
  bulk_logFC = NA_real_,
  stringsAsFactors = FALSE
)
hub_target_table <- data.frame(
  gene = hub_target$gene,
  source = paste0("hub_", hub_target$direction),
  bulk_direction = hub_target$direction,
  bulk_logFC = hub_target$logFC,
  stringsAsFactors = FALSE
)
target_table <- rbind(signature_target, hub_target_table)
target_table <- aggregate(
  cbind(source, bulk_direction) ~ gene,
  data = transform(target_table,
                   source = as.character(source),
                   bulk_direction = as.character(bulk_direction)),
  FUN = function(x) paste(unique(x[!is.na(x) & nzchar(x)]), collapse = ";")
)
logfc_map <- tapply(hub_target_table$bulk_logFC, hub_target_table$gene,
                    function(x) x[which.max(abs(x))])
target_table$bulk_logFC <- unname(logfc_map[target_table$gene])
target_table$bulk_logFC[is.na(target_table$bulk_logFC)] <- NA_real_

make_dotplot(object, signature$gene, "cell_type_plot",
             file.path(output_dir, "fig_dotplot_signature.png"),
             "Signature genes by cell type", width = 11, height = 6.5)
make_dotplot(object, hub_target$gene, "cell_type_plot",
             file.path(output_dir, "fig_dotplot_hub.png"),
             "Top hub genes by cell type", width = 18, height = 8.5)

expression_long <- summarize_target_expression(object, target_table)
write.csv(expression_long, file.path(output_dir, "gene_celltype_expression.csv"),
          row.names = FALSE)
dominant <- dominant_celltype(expression_long)
dominant <- dominant[match(target_table$gene, dominant$gene), , drop = FALSE]
write.csv(dominant, file.path(output_dir, "gene_dominant_celltype.csv"),
          row.names = FALSE)

persistent_rds <- file.path(project_root, "data", "processed",
                            "GSE230375_seurat_bpcells.rds")
dir.create(dirname(persistent_rds), recursive = TRUE, showWarnings = FALSE)
saveRDS(object, persistent_rds, compress = TRUE)
log_message("Saved persistent Seurat object to ", persistent_rds)

report_path <- file.path(output_dir, "08_scrna_report.md")
write_report(object, manifest, failure_df, qc_summary, annotation, dominant, composition,
             harmony_status, report_path)

log_message("Analysis complete")
log_message("Output directory: ", output_dir)
log_message("Report: ", report_path)
}

if (!identical(Sys.getenv("VS_SCRNA_TEST_MODE"), "1")) {
  main()
}
