#!/usr/bin/env Rscript
# 18_figures_v2.R — 重排 + 加料 + 统一风格的投稿图(矢量 PDF + PNG 预览)
# 复用 10_publication_figures.R 的 theme_pub / 配色 / patchwork 思路, 但:
#   (1) 主图重排: Fig1 DEG | Fig2 跨平台 | Fig3 单细胞组成 | Fig4 组成感知分解 | Fig5 GSE216783
#   (2) 给信息量薄的图加料(把相关辅助内容并入主图): Fig1+富集, Fig2+逐基因幅度, Fig3+组成柱
#   (3) 全部数据派生面板用同一 theme_pub 重渲, patchwork 真对齐 + 自动 A/B/C 标签
#   (4) 矢量 cairo_pdf(供 Adobe 调整) + 300dpi PNG(供 docx 内嵌)
# 输出: results/figures_submission/{Fig1..Fig5,FigS1..FigS6}.{pdf,png,tiff}
options(stringsAsFactors = FALSE, width = 120)
suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(grid)
})

root <- Sys.getenv("PROJ_ROOT", getwd())
outdir <- file.path(root, "results", "figures_submission")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

theme_pub <- function(base_size = 8.5) {
  theme_bw(base_size = base_size) +
    theme(text = element_text(family = "sans", color = "black"),
          axis.text = element_text(color = "black", size = 7.5),
          axis.title = element_text(size = 8.5),
          plot.title = element_text(size = 9, face = "bold", hjust = 0),
          legend.title = element_text(size = 8), legend.text = element_text(size = 7),
          strip.text = element_text(size = 7.5, face = "bold"),
          panel.grid.minor = element_blank(), plot.margin = margin(5,5,5,5))
}
theme_set(theme_pub())
module_colors <- c(turquoise = "#1B9E77", yellow = "#D9A400")
grp_colors <- c(Nerve = "#2166AC", VS = "#B2182B")
deg_colors <- c("Downregulated"="#2166AC","Not significant"="#BDBDBD","Upregulated"="#B2182B")
adj_colors <- c("retained"="#B2182B","attenuated"="#9ECAE1","adjusted-only"="#D9A400",
                "unchanged"="#BDBDBD","background"="#8F8F8F")
zh2en <- c("Schwann/肿瘤细胞"="Schwann/tumour","巨噬/髓系"="Myeloid","T细胞"="T cell",
           "成纤维"="Fibroblast","内皮"="Endothelial","周细胞"="Pericyte")
ct_levels <- c("Schwann/tumour","Myeloid","T cell","Fibroblast","Endothelial","Pericyte")
std_cell <- function(x) {
  y <- as.character(x)
  mapped <- unname(zh2en[y])
  ifelse(is.na(mapped), y, mapped)
}
stars <- function(p) ifelse(is.na(p), "NE", ifelse(p<0.001,"***",ifelse(p<0.01,"**",ifelse(p<0.05,"*","ns"))))
rc <- function(p){ f<-file.path(root,p); if(!file.exists(f)) stop("missing ",p); read.csv(f, check.names=FALSE) }
TAG <- function() plot_annotation(tag_levels="A", theme=theme(plot.tag=element_text(size=13,face="bold")))

export_figure <- function(plot, stem, width, height) {
  grDevices::cairo_pdf(file.path(outdir, paste0(stem,".pdf")), width=width, height=height, family="sans")
  print(plot); grDevices::dev.off()
  grDevices::png(file.path(outdir, paste0(stem,".png")), width=width, height=height, units="in",
                 res=300, type="cairo"); print(plot); grDevices::dev.off()
  grDevices::tiff(file.path(outdir, paste0(stem,".tiff")), width=width, height=height, units="in",
                  res=300, compression="lzw", type="cairo"); print(plot); grDevices::dev.off()
  cat(sprintf("  %-7s %.1f x %.1f in\n", stem, width, height))
}
raster_panel <- function(p){  # 将既有 PNG 作为对齐网格中的一个面板(UMAP/dotplot 本就栅格)
  img <- png::readPNG(file.path(root,p)); wrap_elements(full = rasterGrob(img, interpolate=TRUE))
}

cat("== Figure 1: discovery DE (volcano + heatmap + enrichment) ==\n")
deg <- rc("results/deg/deg_full_gene.csv")
deg$category <- factor(ifelse(deg$adj.P.Val<0.05 & deg$logFC>1,"Upregulated",
                       ifelse(deg$adj.P.Val<0.05 & deg$logFC< -1,"Downregulated","Not significant")),
                       levels=names(deg_colors))
deg$ml10 <- -log10(pmax(deg$adj.P.Val,.Machine$double.xmin))
lab <- do.call(rbind, lapply(split(deg[deg$category!="Not significant",], deg$category[deg$category!="Not significant"]),
              function(x) head(x[order(x$adj.P.Val,-abs(x$logFC)),],5)))
p1a <- ggplot(deg, aes(logFC, ml10, color=category)) +
  geom_point(size=0.6, alpha=0.6) +
  geom_vline(xintercept=c(-1,1), linetype="dashed", linewidth=0.3) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", linewidth=0.3) +
  ggrepel::geom_text_repel(data=lab, aes(label=gene), size=2.4, color="black",
                           min.segment.length=0, box.padding=0.25, max.overlaps=Inf, seed=2026) +
  scale_color_manual(values=deg_colors, drop=FALSE) +
  labs(title="Discovery-cohort differential expression",
       x=expression(log[2]~fold~change), y=expression(-log[10]~adjusted~italic(P)), color=NULL) +
  theme(legend.position="bottom")
# heatmap
expr_df <- rc("data/processed/GSE39645_expr_gene.csv"); rownames(expr_df) <- expr_df$gene
expr <- as.matrix(expr_df[, setdiff(names(expr_df),"gene")]); storage.mode(expr) <- "double"
meta <- rc("data/processed/GSE39645_meta.csv"); meta <- meta[match(colnames(expr),meta$sample),]
sig <- deg[deg$adj.P.Val<0.05 & abs(deg$logFC)>1,]
top50 <- intersect(head(sig$gene[order(sig$adj.P.Val,-abs(sig$logFC))],30), rownames(expr))
so <- order(factor(meta$group,levels=c("Control","Tumor")), meta$batch)
z <- t(scale(t(expr[top50, so]))); z[!is.finite(z)] <- 0; z <- pmax(pmin(z,2.5),-2.5)
ann <- data.frame(Group=ifelse(meta$group[so]=="Tumor","VS","Nerve"), Batch=factor(meta$batch[so]))
rownames(ann) <- meta$sample[so]
hm <- pheatmap::pheatmap(z, cluster_rows=TRUE, cluster_cols=FALSE, show_colnames=FALSE,
        annotation_col=ann, annotation_colors=list(Group=grp_colors,
        Batch=setNames(c("#FDD49E","#7A0177")[seq_along(unique(ann$Batch))], levels(ann$Batch))),
        color=colorRampPalette(c("#2166AC","white","#B2182B"))(101), breaks=seq(-2.5,2.5,length.out=102),
        border_color=NA, fontsize=7, fontsize_row=6, main="Top 30 DEGs (row z-score)", silent=TRUE)
p1b <- wrap_elements(full=hm$gtable)
# enrichment (top GO BP up & down)
gu <- rc("results/enrich/GO_BP_up.csv"); gd <- rc("results/enrich/GO_BP_down.csv")
mk <- function(d,dir,n=6){ d<-head(d[order(d$p.adjust),],n); data.frame(Description=d$Description,
        ml10=-log10(d$p.adjust), dir=dir) }
en <- rbind(mk(gu,"Up in VS"), mk(gd,"Down in VS"))
en$Description <- ifelse(nchar(en$Description)>50, paste0(substr(en$Description,1,49),"…"), en$Description)
en$Description <- factor(en$Description, levels=rev(en$Description))
p1c <- ggplot(en, aes(ml10, Description, fill=dir)) +
  geom_col(width=0.7) +
  scale_fill_manual(values=c("Up in VS"="#B2182B","Down in VS"="#2166AC"), name=NULL) +
  labs(title="GO biological-process enrichment", x=expression(-log[10]~adjusted~italic(P)), y=NULL) +
  theme(legend.position="bottom", axis.text.y=element_text(size=6.5))
fig1 <- (p1a | p1b) / p1c + plot_layout(heights=c(1,0.62)) + TAG()
export_figure(fig1, "Fig1", 13.2, 10.0)

cat("== Figure 2: cross-platform reproducibility (scatter + per-gene magnitude) ==\n")
val <- rc("results/validation/hub_logFC_concordance.csv")
val$module <- factor(val$module, levels=c("turquoise","yellow"))
p2a <- ggplot(val, aes(discovery_logFC, validation_logFC, color=module)) +
  geom_hline(yintercept=0, linewidth=0.25, color="grey70") +
  geom_vline(xintercept=0, linewidth=0.25, color="grey70") +
  geom_abline(slope=1, intercept=0, linetype=2, linewidth=0.4, color="grey50") +
  geom_point(size=1.3, alpha=0.7) + geom_smooth(method="lm", se=TRUE, linewidth=0.7, color="black") +
  scale_color_manual(values=module_colors, drop=FALSE) +
  annotate("label", x=-Inf, y=Inf, hjust=-0.05, vjust=1.15, size=2.8, linewidth=0.2,
           label="Pearson r = 0.883\nDirection concordance = 93.4%") +
  labs(title="Cross-platform hub-gene logFC", x=expression(Discovery~log[2]~FC),
       y=expression(Validation~log[2]~FC), color="Module") + coord_cartesian(clip="off")
# per-gene paired magnitude (top |discovery|)
topg <- head(val[order(-abs(val$discovery_logFC)),], 14)
pg <- rbind(data.frame(gene=topg$gene, lfc=topg$discovery_logFC, cohort="Discovery", module=topg$module),
            data.frame(gene=topg$gene, lfc=topg$validation_logFC, cohort="Validation", module=topg$module))
pg$gene <- factor(pg$gene, levels=topg$gene[order(topg$discovery_logFC)])
p2b <- ggplot(pg, aes(lfc, gene, fill=cohort)) +
  geom_col(position=position_dodge(width=0.7), width=0.65) +
  geom_vline(xintercept=0, linewidth=0.25, color="grey60") +
  scale_fill_manual(values=c("Discovery"="#762A83","Validation"="#1B7837"), name=NULL) +
  labs(title="Per-gene logFC in both cohorts", x=expression(log[2]~fold~change), y=NULL) +
  theme(legend.position="bottom", axis.text.y=element_text(size=6.5, face="italic"))
fig2 <- (p2a | p2b) + plot_layout(widths=c(1,0.95)) + TAG()
export_figure(fig2, "Fig2", 12.6, 6.2)

cat("== Figure 3: single-cell composition & localisation (UMAP raster + donor composition) ==\n")
comp <- rc("results/pseudobulk/pseudobulk_sample_celltype_ncells.csv")
comp$within_sample_fraction <- comp$n_cells / ave(comp$n_cells, comp$sample, FUN=sum)
comp$cell_type <- factor(std_cell(comp$cell_type), levels=ct_levels)
comp$group <- if ("group" %in% names(comp)) comp$group else ifelse(startsWith(comp$sample,"N"),"Nerve","VS")
comp$group <- factor(ifelse(comp$group=="Tumor","VS",comp$group), levels=c("Nerve","VS"))
p3d <- ggplot(comp, aes(cell_type, within_sample_fraction, fill=group)) +
  geom_boxplot(position=position_dodge(width=0.75), outlier.size=0.4, width=0.66, linewidth=0.3) +
  geom_point(position=position_jitterdodge(jitter.width=0.15, dodge.width=0.75),
             size=0.6, alpha=0.7) +
  scale_fill_manual(values=grp_colors, name=NULL) +
  labs(title="Cell-type composition (per donor)", x=NULL, y="Within-sample fraction") +
  theme(legend.position="bottom", axis.text.x=element_text(angle=20,hjust=1))
fig3 <- (raster_panel("results/scrna/fig_umap_celltype.png") | raster_panel("results/scrna/fig_umap_group.png")) /
        (raster_panel("results/scrna/fig_dotplot_signature.png") | p3d) +
        plot_layout(heights=c(1,0.95)) + TAG()
export_figure(fig3, "Fig3", 13.0, 10.6)

cat("== Figure 4: composition-aware decomposition (4 re-rendered panels) ==\n")
# A direction concordance (descriptive, no CI)
dc <- rc("results/pseudobulk/dominant_celltype_direction_concordance.csv")
dca <- dc[dc$scope=="by_celltype" & dc$filter_type=="all",]
dca$dominant_celltype_en <- std_cell(dca$dominant_celltype)
dca$lab <- factor(dca$dominant_celltype_en, levels=dca$dominant_celltype_en[order(-dca$direction_concordance)])
ov <- dc$direction_concordance[dc$scope=="overall" & dc$filter_type=="all"]
ovn <- dc$n_direction[dc$scope=="overall" & dc$filter_type=="all"]
p4a <- ggplot(dca, aes(lab, direction_concordance)) +
  geom_col(fill="#2c7fb8", width=0.66) +
  geom_hline(yintercept=ov, linetype=2, color="grey30") +
  geom_text(aes(label=sprintf("%.1f%%\n(n=%d)",100*direction_concordance,n_direction)), vjust=-0.2, size=2.5) +
  annotate("text", x=6.3, y=ov+0.03, hjust=1, size=2.6, label=sprintf("overall %.1f%% (n=%d)",100*ov,ovn)) +
  scale_y_continuous(limits=c(0,1.05), expand=expansion(mult=c(0,.02))) +
  labs(title="Bulk-DEG direction concordance with dominant cell type",
       x=NULL, y="Direction concordance (descriptive)",
       caption="Genes are not independent units; no inferential interval shown") +
  theme(axis.text.x=element_text(angle=20,hjust=1), plot.caption=element_text(size=6.5,color="grey40"))
# B MuSiC fractions GSE39645 + stars
fr <- rc("results/deconv/cell_fractions_GSE39645_music.csv")
frl <- reshape(fr, varying=list(names(fr)[4:9]), v.names="frac", timevar="cell_type",
               times=names(fr)[4:9], direction="long")
frl$cell_type <- factor(std_cell(frl$cell_type), levels=ct_levels)
frl$group <- factor(ifelse(frl$group=="Tumor","VS","Nerve"), levels=c("Nerve","VS"))
gt <- rc("results/deconv/cell_fraction_group_tests.csv")
gt <- gt[gt$cohort=="GSE39645" & gt$method=="MuSiC",]
gt$cell_type <- factor(std_cell(gt$cell_type), levels=ct_levels)
gt$star <- stars(gt$padj)
p4b <- ggplot(frl, aes(cell_type, frac, fill=group)) +
  geom_boxplot(width=0.66, outlier.size=0.4, linewidth=0.3, position=position_dodge(0.72)) +
  geom_text(data=gt, aes(cell_type, 1.02, label=star), inherit.aes=FALSE, size=2.7) +
  scale_fill_manual(values=grp_colors, name=NULL) +
  scale_y_continuous(limits=c(0,1.08)) +
  labs(title="MuSiC relative cell fractions (GSE39645)", x=NULL, y="Estimated relative fraction") +
  theme(axis.text.x=element_text(angle=20,hjust=1), legend.position="bottom")
# C composition-adjusted scatter
ca <- rc("results/deconv/deg_composition_adjusted.csv")
ca$classification[is.na(ca$classification) | ca$classification=="NA"] <- "background"
ca$classification[!(ca$classification %in% names(adj_colors))] <- "background"
ca$classification <- factor(ca$classification, levels=names(adj_colors))
p4c <- ggplot(ca, aes(original_logFC, adjusted_logFC, color=classification)) +
  geom_abline(slope=1, intercept=0, linetype=2, color="grey50") +
  geom_point(size=0.5, alpha=0.5) +
  scale_color_manual(values=adj_colors, name=NULL, drop=FALSE) +
  labs(title="DE before vs after composition adjustment (GSE39645)",
       x=expression(Original~log[2]~FC), y=expression(Adjusted~log[2]~FC)) +
  guides(color=guide_legend(override.aes=list(size=2))) + theme(legend.position="bottom")
# D bulk-DEG classification by composition consistency and within-compartment evidence
ivc <- rc("results/pseudobulk/intrinsic_vs_composition.csv")
class_levels <- c("composition-consistent","both","cell-intrinsic (candidate)","ambiguous/undetermined")
class_colors <- c("composition-consistent"="#5B7C99","both"="#762A83",
                  "cell-intrinsic (candidate)"="#B2182B","ambiguous/undetermined"="#BDBDBD")
class_counts <- as.data.frame(table(factor(ivc$class, levels=class_levels)))
names(class_counts) <- c("class","count")
class_counts$percent <- 100 * class_counts$count / sum(class_counts$count)
class_counts$class <- factor(class_counts$class, levels=rev(class_levels))
p4class <- ggplot(class_counts, aes(count, class, fill=class)) +
  geom_col(width=0.68) +
  geom_text(aes(label=sprintf("%d (%.1f%%)",count,percent)), hjust=-0.08, size=2.6) +
  scale_fill_manual(values=class_colors, guide="none", drop=FALSE) +
  scale_x_continuous(expand=expansion(mult=c(0,0.18))) +
  labs(title="Bulk-DEG composition vs within-compartment classification",
       x="Number of bulk DEGs", y=NULL)
# D representative pseudobulk vs bulk (Schwann/tumour, ALL matched genes — matches legend r=0.179/n=14369)
pbsch <- rc("results/pseudobulk/Schwann_tumor_VS_vs_nerve_DE.csv")[, c("gene","log2FC")]
ivs <- merge(deg[, c("gene","logFC")], pbsch, by="gene")
ivs <- ivs[is.finite(ivs$logFC) & is.finite(ivs$log2FC), ]
names(ivs)[names(ivs)=="logFC"] <- "bulk_logFC"; names(ivs)[names(ivs)=="log2FC"] <- "pb_log2FC_in_dom"
rval <- suppressWarnings(cor(ivs$bulk_logFC, ivs$pb_log2FC_in_dom, use="complete.obs"))
p4d <- ggplot(ivs, aes(bulk_logFC, pb_log2FC_in_dom)) +
  geom_hline(yintercept=0, linewidth=0.2, color="grey75") + geom_vline(xintercept=0, linewidth=0.2, color="grey75") +
  geom_point(size=0.5, alpha=0.4, color="#B2182B") +
  geom_smooth(method="lm", se=TRUE, linewidth=0.7, color="black") +
  annotate("label", x=-Inf, y=Inf, hjust=-0.05, vjust=1.15, size=2.7, linewidth=0.2,
           label=sprintf("Pearson r = %.3f (n=%d)", rval, nrow(ivs))) +
  labs(title="Within-Schwann/tumour pseudobulk vs bulk",
       x=expression(Bulk~log[2]~FC), y=expression(Pseudobulk~log[2]~FC))
fig4 <- (p4a | p4b) / (p4c | p4class) + plot_layout(heights=c(1,1)) + TAG()
export_figure(fig4, "Fig4", 12.8, 10.2)

cat("== Figure 5: GSE216783 reference-sensitivity (3 re-rendered panels) ==\n")
g5a <- rc("results/gse216783/compartment_composition.csv")
g5a$cell_type <- std_cell(g5a$cell_type)
g5a <- g5a[g5a$cell_type %in% ct_levels,]
g5a$cell_type <- factor(g5a$cell_type, levels=g5a$cell_type[order(-g5a$percent)])
p5a <- ggplot(g5a, aes(cell_type, percent, fill=cell_type)) +
  geom_col(width=0.7, show.legend=FALSE) +
  geom_text(aes(label=sprintf("%.1f%%",percent)), vjust=-0.3, size=2.5) +
  scale_fill_brewer(palette="Set2") +
  scale_y_continuous(expand=expansion(mult=c(0,.12))) +
  labs(title="GSE216783 within-VS composition", x=NULL, y="% of cells") +
  theme(axis.text.x=element_text(angle=20,hjust=1))
g5b <- rc("results/gse216783/candidate_compartment_localization.csv")
g5b$dominant_celltype_gse216783 <- std_cell(g5b$dominant_celltype_gse216783)
g5b$gene <- factor(g5b$gene, levels=g5b$gene)
p5b <- ggplot(g5b, aes(gene, 1, fill=concordant)) +
  geom_tile(color="white", linewidth=1) +
  geom_text(aes(label=dominant_celltype_gse216783), angle=90, size=2.2, color="grey15", fontface="bold") +
  scale_fill_manual(values=c("TRUE"="#74A9A1","FALSE"="#C9A9A6"), name="Concordant\nwith GSE230375") +
  labs(title="Candidate compartment localisation (9/11 concordant)", x=NULL, y=NULL) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.text.x=element_text(angle=45,hjust=1,face="italic"), panel.grid=element_blank())
g5c <- rc("results/gse216783/reference_comparison_deconv.csv")
g5c$cell_type <- factor(std_cell(g5c$cell_type), levels=ct_levels)
p5c <- ggplot(g5c, aes(cell_type, spearman_vs_GSE230375ref, fill=cohort)) +
  geom_col(data=g5c[is.finite(g5c$spearman_vs_GSE230375ref),],
           position=position_dodge(0.72), width=0.65) +
  geom_text(data=g5c[!is.finite(g5c$spearman_vs_GSE230375ref),],
            aes(x=cell_type, y=0.04, label="NE", group=cohort),
            position=position_dodge(0.72), size=2.3, inherit.aes=FALSE) +
  geom_hline(yintercept=0.7, linetype=2, color="grey40") +
  scale_fill_manual(values=c("GSE39645"="#4575B4","GSE108524"="#D9A400"), name=NULL) +
  scale_y_continuous(limits=c(0,1)) +
  labs(title="2nd-reference vs GSE230375-reference deconvolution",
       x=NULL, y="Spearman r (sample-level)") +
  theme(axis.text.x=element_text(angle=20,hjust=1), legend.position="bottom")
fig5 <- (p5a | p5c) / p5b + plot_layout(heights=c(1,0.7)) + TAG()
export_figure(fig5, "Fig5", 12.6, 8.4)

cat("== Supplementary S1 (WGCNA), S2 (immune) ==\n")
mt <- rc("results/wgcna/module_trait_cor.csv")
mt$module <- factor(mt$module, levels=rev(mt$module[order(mt$cor)]))
mt$label <- sprintf("r = %.2f\n%s", mt$cor, format.pval(mt$p, digits=2, eps=0.001))
ps1a <- ggplot(mt, aes("VS", module, fill=cor)) + geom_tile(color="white", linewidth=0.5) +
  geom_text(aes(label=label), size=2.6) +
  scale_fill_steps2(low="#2166AC", mid="white", high="#B2182B", midpoint=0, limits=c(-1,1), n.breaks=7, name="Correlation") +
  labs(title="Module–trait association", x=NULL, y="Module") +
  theme(panel.grid=element_blank(), axis.text.x=element_text(face="bold"))
me <- rc("results/hub/module_gene_metrics.csv"); me <- me[me$module %in% c("turquoise","yellow"),]
me$module <- factor(me$module, levels=c("turquoise","yellow"))
mk_kme <- function(d, mod) {
  ggplot(d[d$module==mod,], aes(kME, GS)) +
    geom_point(size=0.7, alpha=0.55, color=module_colors[[mod]]) +
    geom_smooth(method="lm", se=TRUE, linewidth=0.6, color="black") +
    labs(title=sprintf("%s module: kME vs GS", mod),
         x="Module membership (kME)", y="Gene significance (GS)") +
    theme(legend.position="none")
}
ps1b <- mk_kme(me, "turquoise")
ps1c <- mk_kme(me, "yellow") +
  theme(legend.position="none")
figs1 <- (ps1a | ps1b | ps1c) + plot_layout(widths=c(0.8,1,1)) + TAG()
export_figure(figs1, "FigS1", 14.2, 5.6)

scores <- rc("results/immune/ssgsea_scores.csv"); idiff <- rc("results/immune/immune_diff_VS_vs_control.csv")
icor <- rc("results/immune/signature_immune_cor.csv")
sc <- idiff$cell_type[idiff$p_adj<0.05]
sl <- reshape2::melt(scores[scores$cell_type %in% sc,], id.vars="cell_type", variable.name="sample", value.name="score")
sl <- merge(sl, meta[,c("sample","group")], by="sample")
sl$group <- factor(ifelse(sl$group=="Tumor","VS","Nerve"), levels=c("Nerve","VS"))
fl <- setNames(paste0(idiff$cell_type[match(sc,idiff$cell_type)]," ",stars(idiff$p_adj[match(sc,idiff$cell_type)])), sc)
sl$cell_type <- factor(sl$cell_type, levels=sc)
ps2a <- ggplot(sl, aes(group, score, fill=group)) +
  geom_boxplot(width=0.6, outlier.shape=NA, linewidth=0.3) + geom_jitter(width=0.13, size=0.5, alpha=0.6) +
  facet_wrap(~cell_type, scales="free_y", ncol=3, labeller=as_labeller(fl)) +
  scale_fill_manual(values=grp_colors) +
  labs(title="Differential ssGSEA immune/stromal scores", x=NULL, y="ssGSEA score", fill=NULL) +
  theme(legend.position="bottom", axis.text.x=element_text(angle=30,hjust=1))
rho <- reshape2::dcast(icor, gene~cell_type, value.var="spearman_rho")
pad <- reshape2::dcast(icor, gene~cell_type, value.var="p_adj")
rm_ <- as.matrix(rho[,-1]); rownames(rm_) <- rho$gene
pm_ <- as.matrix(pad[,-1]); rownames(pm_) <- pad$gene; pm_ <- pm_[rownames(rm_),colnames(rm_)]
hm5 <- pheatmap::pheatmap(rm_, cluster_rows=TRUE, cluster_cols=TRUE,
        color=colorRampPalette(c("#2166AC","white","#B2182B"))(101), breaks=seq(-1,1,length.out=102),
        border_color="white", display_numbers=matrix(stars(pm_),nrow=nrow(pm_),dimnames=dimnames(pm_)),
        number_color="black", fontsize_number=6, fontsize=7, fontsize_row=7, fontsize_col=6,
        angle_col=45, main="Transcript–immune score correlations", silent=TRUE)
ps2b <- wrap_elements(full=hm5$gtable)
figs2 <- (ps2a | ps2b) + plot_layout(widths=c(1.35,1)) + TAG()
export_figure(figs2, "FigS2", 14.0, 8.4)

cat("== Supplementary S3 (pseudobulk volcanoes), S4 (bulk/pseudobulk concordance), S5 (validation deconv), S6 (module preservation) ==\n")
volc <- function(file, title){
  d <- rc(file)
  d <- d[is.finite(d$log2FC) & is.finite(d$padj),]
  d$cat <- ifelse(d$padj<0.05 & d$log2FC>1,"Up",ifelse(d$padj<0.05 & d$log2FC< -1,"Down","ns"))
  d$ml10 <- -log10(pmax(d$padj,1e-300))
  ggplot(d, aes(log2FC, ml10, color=cat)) + geom_point(size=0.4, alpha=0.5) +
    geom_vline(xintercept=c(-1,1), linetype=2, linewidth=0.25) + geom_hline(yintercept=-log10(0.05), linetype=2, linewidth=0.25) +
    scale_color_manual(values=c("Up"="#B2182B","Down"="#2166AC","ns"="#BDBDBD"), guide="none") +
    labs(title=title, x=expression(log[2]~FC), y=expression(-log[10]~adjP)) + theme(plot.title=element_text(size=8))
}
pv <- list(
  volc("results/pseudobulk/Schwann_tumor_VS_vs_nerve_DE.csv","Schwann/tumour"),
  volc("results/pseudobulk/Myeloid_VS_vs_nerve_DE.csv","Myeloid"),
  volc("results/pseudobulk/T_cell_VS_vs_nerve_DE.csv","T cell"),
  volc("results/pseudobulk/Fibroblast_VS_vs_nerve_DE.csv","Fibroblast"),
  volc("results/pseudobulk/Endothelial_VS_vs_nerve_DE.csv","Endothelial"),
  volc("results/pseudobulk/Pericyte_VS_vs_nerve_DE.csv","Pericyte"))
figs3 <- wrap_plots(pv, ncol=3) + TAG()
export_figure(figs3, "FigS3", 12.0, 7.2)

fr8 <- rc("results/deconv/cell_fractions_GSE108524_music.csv")
fr8l <- reshape(fr8, varying=list(names(fr8)[4:9]), v.names="frac", timevar="cell_type",
                times=names(fr8)[4:9], direction="long")
fr8l$cell_type <- factor(std_cell(fr8l$cell_type), levels=ct_levels)
fr8l$group <- factor(ifelse(fr8l$group=="Tumor","VS","Nerve"), levels=c("Nerve","VS"))
gt8 <- rc("results/deconv/cell_fraction_group_tests.csv"); gt8 <- gt8[gt8$cohort=="GSE108524" & gt8$method=="MuSiC",]
gt8$cell_type <- factor(std_cell(gt8$cell_type), levels=ct_levels); gt8$star <- stars(gt8$padj)
ps6 <- ggplot(fr8l, aes(cell_type, frac, fill=group)) +
  geom_boxplot(width=0.66, outlier.size=0.4, linewidth=0.3, position=position_dodge(0.72)) +
  geom_text(data=gt8, aes(cell_type, 1.02, label=star), inherit.aes=FALSE, size=2.7) +
  scale_fill_manual(values=grp_colors, name=NULL) + scale_y_continuous(limits=c(0,1.08)) +
  labs(title="MuSiC relative cell fractions (GSE108524, 4 nerve controls)", x=NULL, y="Estimated relative fraction") +
  theme(axis.text.x=element_text(angle=20,hjust=1), legend.position="bottom")
figs5 <- ps6 + TAG()
export_figure(figs5, "FigS5", 8.0, 5.2)

cat("== Supplementary S4 (bulk vs pseudobulk concordance scatters, raster + representative scatter) ==\n")
figs4 <- (raster_panel("results/pseudobulk/fig_corr_Schwann_tumor_vs_bulk.png") |
          raster_panel("results/pseudobulk/fig_corr_Myeloid_vs_bulk.png") |
          raster_panel("results/pseudobulk/fig_corr_Fibroblast_vs_bulk.png")) /
         (raster_panel("results/pseudobulk/fig_corr_sensitivity_Schwann_tumor.png") |
          raster_panel("results/pseudobulk/fig_corr_sensitivity_Myeloid.png") |
          raster_panel("results/pseudobulk/fig_corr_sensitivity_Fibroblast.png")) /
         (plot_spacer() | p4d | plot_spacer()) + plot_layout(heights=c(1,1,1)) + TAG()
export_figure(figs4, "FigS4", 13.5, 11.5)

cat("== Supplementary S6 (WGCNA module preservation) ==\n")
mp <- rc("results/wgcna/module_preservation.csv")
mp$module <- factor(mp$module, levels=mp$module[order(mp$moduleSize)])
mp$module_class <- ifelse(as.character(mp$module) %in% names(module_colors), as.character(mp$module), "other")
mp_cols <- c(module_colors, other="#9E9E9E")
ps6a <- ggplot(mp, aes(moduleSize, medianRank, label=module, color=module_class)) +
  geom_point(size=2.2) +
  ggrepel::geom_text_repel(size=2.4, min.segment.length=0, max.overlaps=Inf, seed=2026) +
  scale_x_log10() +
  scale_color_manual(values=mp_cols, guide="none") +
  labs(title="Module preservation median rank",
       x="Module size (log scale)", y="Median rank (lower is better)")
ps6b <- ggplot(mp, aes(moduleSize, Zsummary, label=module, color=module_class)) +
  geom_hline(yintercept=c(2,10), linetype="dashed", linewidth=0.35, color=c("#2166AC","#1B7837")) +
  geom_point(size=2.2) +
  ggrepel::geom_text_repel(size=2.4, min.segment.length=0, max.overlaps=Inf, seed=2026) +
  scale_x_log10() +
  scale_color_manual(values=mp_cols, guide="none") +
  labs(title="Module preservation Zsummary",
       x="Module size (log scale)", y="Zsummary")
figs6 <- (ps6a | ps6b) + TAG()
export_figure(figs6, "FigS6", 11.0, 5.4)

cat("DONE ->", outdir, "\n")
