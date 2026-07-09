#!/usr/bin/env Rscript
# 02_annotate_deg_enrich.R — probe->gene 注释 + 基因级 DEG(batch-adj) + 富集
suppressMessages({
  library(data.table); library(limma); library(ggplot2); library(ggrepel); library(pheatmap)
  library(hugene10sttranscriptcluster.db); library(clusterProfiler); library(org.Hs.eg.db)
})
set.seed(42)
root <- Sys.getenv("PROJ_ROOT", getwd())
pd <- file.path(root,"data/processed"); dir.create(pd,showWarnings=FALSE,recursive=TRUE)
dd <- file.path(root,"results/deg"); dir.create(dd,showWarnings=FALSE,recursive=TRUE)
ed <- file.path(root,"results/enrich"); dir.create(ed,showWarnings=FALSE,recursive=TRUE)
mf <- file.path(root,"data/raw/GSE39645_matrix.txt.gz")

## --- parse series matrix (group, batch, expr) ---
L <- readLines(gzfile(mf)); gl <- function(t) L[grepl(paste0("^",t),L)]
samp_titles <- gsub('"','',strsplit(gl("!Sample_geo_accession"),"\t")[[1]][-1])
chars <- lapply(gl("!Sample_characteristics_ch1"), function(x) gsub('"','',strsplit(x,"\t")[[1]][-1]))
group <- factor(sub(".*group:\\s*","", chars[[ which(sapply(chars,function(x)any(grepl("group:",x))))[1] ]]))
bidx <- which(sapply(chars,function(x)any(grepl("batch:",x))))
batch <- factor(sub(".*batch:\\s*","", chars[[bidx[1]]]))
beg <- grep("^!series_matrix_table_begin",L); end <- grep("^!series_matrix_table_end",L)
dt <- fread(text=paste(L[(beg+1):(end-1)],collapse="\n"),header=TRUE); setnames(dt,1,"probe")
expr <- as.matrix(dt[,-1]); rownames(expr) <- as.character(dt$probe)
colnames(expr) <- gsub('"','',colnames(expr))
if(max(expr,na.rm=TRUE)>50) expr <- log2(expr+1)
expr <- expr[complete.cases(expr),]
cat(sprintf("probes=%d samples=%d | Control=%d Tumor=%d\n",nrow(expr),ncol(expr),sum(group=="Control"),sum(group=="Tumor")))

## --- probe -> symbol, collapse by max mean ---
map <- AnnotationDbi::select(hugene10sttranscriptcluster.db, keys=rownames(expr),
        columns=c("SYMBOL"), keytype="PROBEID")
map <- map[!is.na(map$SYMBOL) & !duplicated(map$PROBEID),]
expr2 <- expr[map$PROBEID,]; rownames(expr2) <- map$SYMBOL
# 多 probe→同 symbol: 取平均表达最高者
ord <- order(rowMeans(expr2), decreasing=TRUE)
expr2 <- expr2[ord,]; expr2 <- expr2[!duplicated(rownames(expr2)),]
cat(sprintf("genes after collapse: %d\n", nrow(expr2)))
fwrite(data.table(gene=rownames(expr2), expr2), file.path(pd,"GSE39645_expr_gene.csv"))
fwrite(data.table(sample=colnames(expr2), group=as.character(group), batch=as.character(batch)),
       file.path(pd,"GSE39645_meta.csv"))

## --- gene-level DEG (batch-adjusted) ---
group <- relevel(group, ref="Control")
design <- model.matrix(~0+group+batch); colnames(design)[1:2] <- c("Control","Tumor")
fit <- eBayes(contrasts.fit(lmFit(expr2,design), makeContrasts(Tumor-Control,levels=design)),
              trend=TRUE, robust=TRUE)
tt <- topTable(fit, number=Inf, sort.by="P"); tt$gene <- rownames(tt)
fwrite(tt[,c("gene","logFC","AveExpr","t","P.Value","adj.P.Val","B")], file.path(dd,"deg_full_gene.csv"))
sig <- subset(tt, adj.P.Val<0.05 & abs(logFC)>1)
fwrite(sig[,c("gene","logFC","P.Value","adj.P.Val")], file.path(dd,"deg_sig_fdr05_fc1.csv"))
cat(sprintf("gene DEG FDR<0.05&|logFC|>1: %d (up %d / down %d)\n", nrow(sig), sum(sig$logFC>0), sum(sig$logFC<0)))

## --- volcano ---
tt$sig <- ifelse(tt$adj.P.Val<0.05 & abs(tt$logFC)>1, ifelse(tt$logFC>0,"Up","Down"),"NS")
lab <- head(tt[order(tt$P.Value),],15)
ggsave(file.path(dd,"fig_volcano.png"),
  ggplot(tt,aes(logFC,-log10(P.Value),color=sig))+geom_point(alpha=.5,size=1)+
    scale_color_manual(values=c(Up="#d62728",Down="#1f77b4",NS="grey80"))+
    geom_text_repel(data=lab,aes(label=gene),size=3,color="black",max.overlaps=20)+
    geom_vline(xintercept=c(-1,1),lty=2,color="grey50")+geom_hline(yintercept=-log10(0.05),lty=2,color="grey50")+
    labs(title="VS (Tumor) vs Control nerve — GSE39645", x="log2FC", y="-log10(P)")+theme_bw(),
  width=7,height=5.5,dpi=150)

## --- heatmap top50 by |logFC| among sig ---
hg <- head(sig$gene[order(-abs(sig$logFC))], 50)
ann <- data.frame(group=group,batch=batch); rownames(ann)<-colnames(expr2)
pheatmap(t(scale(t(expr2[hg,]))), annotation_col=ann, show_colnames=FALSE, fontsize_row=6,
         filename=file.path(dd,"fig_heatmap_top50.png"), width=8, height=9)

## --- enrichment (ORA: up & down separately) ---
to_entrez <- function(g) bitr(g, "SYMBOL","ENTREZID", org.Hs.eg.db)$ENTREZID
run_ora <- function(genes, tag){
  e <- tryCatch(to_entrez(genes), error=function(x)NULL); if(is.null(e)||!length(e)) return(invisible())
  go <- enrichGO(e, org.Hs.eg.db, ont="BP", pvalueCutoff=0.05, qvalueCutoff=0.1, readable=TRUE)
  kg <- tryCatch(enrichKEGG(e, pvalueCutoff=0.1), error=function(x)NULL)
  if(!is.null(go)&&nrow(as.data.frame(go))) fwrite(as.data.frame(go), file.path(ed,paste0("GO_BP_",tag,".csv")))
  if(!is.null(kg)&&nrow(as.data.frame(kg))) fwrite(as.data.frame(kg), file.path(ed,paste0("KEGG_",tag,".csv")))
  cat(sprintf("ORA %s: GO_BP=%d KEGG=%d\n", tag,
      ifelse(is.null(go),0,nrow(as.data.frame(go))), ifelse(is.null(kg),0,nrow(as.data.frame(kg)))))
}
run_ora(subset(sig,logFC>0)$gene, "up")
run_ora(subset(sig,logFC<0)$gene, "down")

## --- GSEA (GO:BP) on t-stat ranking ---
ranks <- tt$t; names(ranks) <- mapIds(org.Hs.eg.db, tt$gene, "ENTREZID","SYMBOL")
ranks <- ranks[!is.na(names(ranks))]; ranks <- sort(ranks, decreasing=TRUE)
gse <- tryCatch(gseGO(ranks, OrgDb=org.Hs.eg.db, ont="BP", pvalueCutoff=0.05, verbose=FALSE, eps=0), error=function(x)NULL)
if(!is.null(gse)&&nrow(as.data.frame(gse))){
  g <- as.data.frame(gse); fwrite(g, file.path(ed,"GSEA_GO_BP.csv"))
  cat(sprintf("GSEA GO:BP significant: %d ; top: %s\n", nrow(g), paste(head(g$Description,5),collapse=" | ")))
}
cat("DONE 02\n")
