#!/usr/bin/env Rscript
# 15_candidate_matrix.R
# 跨证据候选映射表: 把 bulk DEG / 跨队列复现 / WGCNA hub / pseudobulk 细胞类型内候选证据 /
#   组成校正保留 五层描述性证据 join 到候选基因。用于 compartment-aware cross-evidence mapping。
# 措辞: in-silico 描述性证据汇总, 不构成机制、因果、标志物或客观排名。每列标注来源。

options(stringsAsFactors = FALSE, warn = 1)
root <- Sys.getenv("PROJ_ROOT", getwd())
outdir <- file.path(root, "results")
rd <- function(p) read.csv(file.path(root, p), check.names = FALSE)

deg  <- rd("results/deg/deg_full_gene.csv")[, c("gene","logFC","adj.P.Val")]
names(deg) <- c("gene","bulk_logFC","bulk_adjP")
hub  <- rd("results/hub/hub_all.csv")[, c("gene","module","direction","kME")]
val  <- rd("results/validation/hub_logFC_concordance.csv")[, c("gene","validation_logFC","direction_concordant")]
ivc  <- rd("results/pseudobulk/intrinsic_vs_composition.csv")[, c("gene","dominant_celltype","pb_log2FC_in_dom","pb_p_in_dom","class")]
names(ivc) <- c("gene","pb_dominant_celltype","pb_log2FC","pb_p","pb_class")
adj  <- rd("results/deconv/deg_composition_adjusted.csv")[, c("gene","classification")]
names(adj)[2] <- "composition_adj_class"

# 全基因跨队列方向一致 (GSE108524 logFC = mean(Tumor)-mean(Control); 覆盖非 hub 候选)
v_expr <- rd("data/processed/GSE108524_expr_gene.csv")
v_meta <- rd("data/processed/GSE108524_meta.csv")
vm <- as.matrix(v_expr[, -1]); rownames(vm) <- v_expr[[1]]
tum <- v_meta$sample[v_meta$group == "Tumor"]; ctl <- v_meta$sample[v_meta$group == "Control"]
tum <- intersect(tum, colnames(vm)); ctl <- intersect(ctl, colnames(vm))
val_full <- data.frame(gene = rownames(vm),
                       validation_logFC_full = rowMeans(vm[, tum, drop=FALSE]) - rowMeans(vm[, ctl, drop=FALSE]))

# 主导细胞中文->英文
ct_map <- c("Schwann/肿瘤细胞"="Schwann/tumour","巨噬/髓系"="Myeloid","T细胞"="T cell",
            "成纤维"="Fibroblast","内皮"="Endothelial","周细胞"="Pericyte")
ivc$pb_dominant_celltype <- ifelse(ivc$pb_dominant_celltype %in% names(ct_map),
                                   ct_map[ivc$pb_dominant_celltype], ivc$pb_dominant_celltype)

# 全表 join:
# - 以 DEG-constrained WGCNA hubs 为主体资源；
# - 另加入文中预先声明的 illustrative compartment transcripts，确保 NRG1/L1CAM/FCGBP/PRRX1
#   等跨区室示例可在同一 descriptive map 中展示；
# - 不再使用 LASSO classifier 作为候选池来源，避免 feature-selection leakage 成为投稿攻击面。
illustrative <- c("NRG1","L1CAM","NCAM2","FCGBP","RRM2B","SLC16A7",
                  "PRRX1","HRH1","AR","SHOX2","PDGFRA")
cand <- sort(unique(c(hub$gene, illustrative)))
m <- data.frame(gene = cand)
for (tb in list(deg, hub, val, ivc, adj, val_full)) m <- merge(m, tb, by = "gene", all.x = TRUE)

# 全基因跨队列方向一致 (优先用全基因计算, 覆盖非 hub 候选)
m$cross_cohort_concordant <- !is.na(m$validation_logFC_full) & !is.na(m$bulk_logFC) &
  sign(m$validation_logFC_full) == sign(m$bulk_logFC)

# 描述性证据层计数 (透明, 仅作表格浏览辅助, 非统计权重/排名)
m$ev_bulk_sig    <- !is.na(m$bulk_adjP) & m$bulk_adjP < 0.05 & abs(m$bulk_logFC) > 1
m$ev_cross_conc  <- m$cross_cohort_concordant
m$ev_hub         <- !is.na(m$module)
m$ev_pb_intrinsic<- m$pb_class %in% c("both","cell-intrinsic (candidate)")
m$ev_adj_retain  <- m$composition_adj_class %in% "retained"
m$evidence_layer_count <- rowSums(m[, c("ev_bulk_sig","ev_cross_conc","ev_hub","ev_pb_intrinsic","ev_adj_retain")])
m <- m[order(-m$evidence_layer_count, m$bulk_adjP), ]
write.csv(m, file.path(outdir, "compartment_cross_evidence_map.csv"), row.names = FALSE)
write.csv(m, file.path(outdir, "candidate_priority_matrix.csv"), row.names = FALSE)
cat("Full cross-evidence map rows:", nrow(m), " -> results/compartment_cross_evidence_map.csv\n")
cat("Legacy compatibility copy -> results/candidate_priority_matrix.csv\n")

# illustrative transcripts 聚焦表 (markdown, 供手稿/补充表说明)
f <- m[m$gene %in% illustrative, ]
f <- f[order(-f$evidence_layer_count, f$bulk_adjP), ]
col <- c("gene","bulk_logFC","bulk_adjP","validation_logFC_full","cross_cohort_concordant",
         "module","pb_dominant_celltype","pb_log2FC","pb_p","pb_class",
         "composition_adj_class","evidence_layer_count")
f2 <- f[, col]
num <- sapply(f2, is.numeric); f2[num] <- lapply(f2[num], function(z) signif(z, 3))
md <- c(
  "# Compartment-aware cross-evidence map (descriptive; not an additive ranking)",
  "",
  paste0("> Supplementary Table S1. 下表为文中 illustrative transcripts 的子集; 完整 map 见 `results/compartment_cross_evidence_map.csv` (", nrow(m), " candidates = 346 DEG-constrained WGCNA hubs plus pre-specified illustrative compartment transcripts)。"),
  paste0("> 读法限制: `evidence_layer_count` 是满足层数(0-5) 的描述性计数, 非加权排名/统计权重/独立证据计数; 证据层不独立(hub 含 DEG; pseudobulk 与 deconvolution 同源 GSE230375); ", sum(m$evidence_layer_count==5), " 基因满 5 层、", sum(m$evidence_layer_count==4), " 基因满 4 层。NRG1/L1CAM/FCGBP/PRRX1 仅作跨区室 illustrative 示例, 非唯一/客观筛选。in-silico 描述性汇总, 不构成机制/因果/标志物。"),
  "",
  paste0("| ", paste(names(f2), collapse=" | "), " |"),
  paste0("|", paste(rep("---", ncol(f2)), collapse="|"), "|"),
  apply(f2, 1, function(z) paste0("| ", paste(z, collapse=" | "), " |")),
  "",
  "## 证据层来源 (绝对路径)",
  "- bulk_logFC/bulk_adjP: `./results/deg/deg_full_gene.csv` (limma, batch-adj)",
  "- validation_logFC_full/cross_cohort_concordant: GSE108524 全基因 logFC=mean(Tumor)-mean(Control) (`./data/processed/GSE108524_expr_gene.csv`); 方向与发现集一致性, 覆盖非 hub 候选",
  "- module/kME: `./results/hub/hub_all.csv` (WGCNA)",
  "- pb_dominant_celltype/pb_log2FC/pb_p/pb_class: `./results/pseudobulk/intrinsic_vs_composition.csv` (GSE230375 pseudobulk)",
  "- composition_adj_class: `./results/deconv/deg_composition_adjusted.csv` (保守组成校正; retained=校正后仍达阈值)",
  "",
  "## 说明",
  "- `evidence_layer_count` 高仅表示多层 in-silico 证据方向一致, 不等于机制证实; pb_class/校正均为保守判据(nerve n=2; 校正可能误伤真信号)。",
  paste0("- `evidence_layer_count` 不是排序依据: ", sum(m$evidence_layer_count==5), " 基因满 5 层、", sum(m$evidence_layer_count==4), " 基因满 4 层; 表中 NRG1/L1CAM/FCGBP/PRRX1 为 illustrative 跨区室示例, 不主张唯一性或客观最高优先级。"),
  "- LASSO classifier 结果不再作为候选池来源；既有 signature 输出仅作为分析归档，不进入主文推断。",
  "- nominal pb_p<0.05 为未校正筛查(无多重校正; nerve n=2), 仅作描述性证据层, 不作统计权重。"
)
writeLines(md, file.path(root, "docs", "candidate_priority_matrix.md"), useBytes = TRUE)
cat("Illustrative transcript table -> docs/candidate_priority_matrix.md\n")
print(f2[, c("gene","pb_dominant_celltype","pb_class","composition_adj_class","cross_cohort_concordant","evidence_layer_count")])
