#!/usr/bin/env Rscript
# 04_wgcna_hub.R — WGCNA 模块↔肿瘤 + hub 基因(kME ∩ GS ∩ DEG)
suppressMessages({library(data.table); library(limma); library(WGCNA); library(pheatmap)})
options(stringsAsFactors=FALSE); disableWGCNAThreads(); cor <- WGCNA::cor; set.seed(42)
root <- Sys.getenv("PROJ_ROOT", getwd())
wd <- file.path(root,"results/wgcna"); dir.create(wd,showWarnings=FALSE,recursive=TRUE)

expr <- as.data.frame(fread(file.path(root,"data/processed/GSE39645_expr_gene.csv")))
rownames(expr) <- expr$gene; expr$gene <- NULL; expr <- as.matrix(expr)
meta <- fread(file.path(root,"data/processed/GSE39645_meta.csv"))
meta <- meta[match(colnames(expr), meta$sample)]
group <- factor(meta$group, levels=c("Control","Tumor")); batch <- factor(meta$batch)

## 去批次(保留 group)用于网络与可视化
expradj <- removeBatchEffect(expr, batch=batch, design=model.matrix(~group))

## top MAD genes
mad <- apply(expradj,1,mad); sel <- names(sort(mad,decreasing=TRUE))[1:min(5000,sum(mad>0))]
datExpr <- t(expradj[sel,])
gsg <- goodSamplesGenes(datExpr, verbose=0); datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
cat(sprintf("WGCNA input: %d samples x %d genes\n", nrow(datExpr), ncol(datExpr)))

## soft threshold
powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector=powers, networkType="signed", verbose=0)
pw <- sft$powerEstimate; if(is.na(pw)) pw <- 12
cat(sprintf("soft power = %d (R2=%.2f)\n", pw, sft$fitIndices$SFT.R.sq[sft$fitIndices$Power==pw]))

net <- blockwiseModules(datExpr, power=pw, networkType="signed", TOMType="signed",
        minModuleSize=30, mergeCutHeight=0.25, numericLabels=TRUE, maxBlockSize=6000,
        saveTOMs=FALSE, verbose=0)
moduleColors <- labels2colors(net$colors)
cat("modules:\n"); print(table(moduleColors))

## module-trait correlation
MEs <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
trait <- as.numeric(group[gsg$goodSamples]=="Tumor")
mtCor <- cor(MEs, trait, use="p"); mtP <- corPvalueStudent(mtCor, nrow(datExpr))
mt <- data.table(module=sub("^ME","",rownames(mtCor)), cor=as.numeric(mtCor), p=as.numeric(mtP))
mt <- mt[order(-abs(cor))]; fwrite(mt, file.path(wd,"module_trait_cor.csv"))
cat("top module-trait:\n"); print(head(mt,6))

## heatmap
png(file.path(wd,"fig_module_trait.png"), width=600, height=900, res=130)
labeledHeatmap(Matrix=mtCor, xLabels="Tumor", yLabels=rownames(mtCor),
  ySymbols=rownames(mtCor), colorLabels=FALSE, colors=blueWhiteRed(50),
  textMatrix=paste0(signif(mtCor,2),"\n(",signif(mtP,1),")"), setStdMargins=FALSE,
  cex.text=0.6, zlim=c(-1,1), main="Module-Tumor association")
dev.off()

## hub genes in top tumor-associated module
topmod <- mt$module[1]
cat(sprintf("top tumor module: %s (cor=%.2f, p=%.1e)\n", topmod, mt$cor[1], mt$p[1]))
kME <- signedKME(datExpr, MEs); colnames(kME) <- sub("^kME","",colnames(kME))
GS <- as.numeric(cor(datExpr, trait, use="p")); names(GS) <- colnames(datExpr)
inmod <- names(moduleColors)[moduleColors==topmod] # note: moduleColors named by gene
inmod <- colnames(datExpr)[moduleColors==topmod]
hub_tab <- data.table(gene=inmod, kME=kME[inmod, topmod], GS=GS[inmod])
# 叠加 DEG
deg <- fread(file.path(root,"results/deg/deg_full_gene.csv"))
hub_tab <- merge(hub_tab, deg[,.(gene,logFC,adj.P.Val)], by="gene", all.x=TRUE)
hub_tab <- hub_tab[order(-abs(kME))]
fwrite(hub_tab, file.path(wd,"module_genes_topmod.csv"))
hub <- hub_tab[abs(kME)>0.8 & abs(GS)>0.5 & adj.P.Val<0.05 & abs(logFC)>1]
fwrite(hub, file.path(wd,"hub_candidates.csv"))
cat(sprintf("hub candidates (|kME|>0.8 & |GS|>0.5 & DEG): %d\n", nrow(hub)))
cat("top hubs:\n"); print(head(hub[order(-abs(kME))][,.(gene,kME,GS,logFC,adj.P.Val)],20))
saveRDS(list(net=net,moduleColors=moduleColors,datExpr=datExpr,trait=trait,topmod=topmod), file.path(wd,"wgcna_obj.rds"))
cat("DONE 04\n")
