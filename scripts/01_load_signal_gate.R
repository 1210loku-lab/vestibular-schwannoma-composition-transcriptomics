#!/usr/bin/env Rscript
# 01_load_signal_gate.R — GSE39645 (VS vs control nerve) 载入 + 信号门(DEG+置换检验)
# 目的: 在投入完整管线前，客观确认数据有真实组别信号(吸取 GSE186505 教训)
suppressMessages({library(data.table); library(limma)})
set.seed(42)
root <- Sys.getenv("PROJ_ROOT", getwd())
outdir <- file.path(root,"results/signal"); dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
mf <- file.path(root,"data/raw/GSE39645_matrix.txt.gz")

## --- parse series matrix ---
L <- readLines(gzfile(mf))
get_line <- function(tag) L[grepl(paste0("^",tag), L)]
titles <- strsplit(get_line("!Sample_title"),"\t")[[1]][-1]
titles <- gsub('"','',titles)
char_lines <- get_line("!Sample_characteristics_ch1")
parse_char <- function(line){ v <- gsub('"','',strsplit(line,"\t")[[1]][-1]); v }
chars <- lapply(char_lines, parse_char)
# 找 group 与 batch 行
grp_row  <- chars[[ which(sapply(chars, function(x) any(grepl("group:", x)))) [1] ]]
bat_idx  <- which(sapply(chars, function(x) any(grepl("batch:", x))))
group <- factor(sub(".*group:\\s*","", grp_row))
batch <- if(length(bat_idx)) factor(sub(".*batch:\\s*","", chars[[bat_idx[1]]])) else NULL
cat("group table:\n"); print(table(group))
if(!is.null(batch)){cat("batch table:\n"); print(table(batch)); cat("group x batch:\n"); print(table(group,batch))}

## --- expression table ---
beg <- grep("^!series_matrix_table_begin", L); end <- grep("^!series_matrix_table_end", L)
dt <- fread(text=paste(L[(beg+1):(end-1)], collapse="\n"), header=TRUE)
setnames(dt, 1, "ID_REF")
expr <- as.matrix(dt[,-1]); rownames(expr) <- dt$ID_REF
colnames(expr) <- gsub('"','',colnames(expr))
cat(sprintf("expr matrix: %d probes x %d samples\n", nrow(expr), ncol(expr)))
stopifnot(length(group)==ncol(expr))
# 值域判断是否已 log
rng <- range(expr, na.rm=TRUE); cat(sprintf("value range: %.2f .. %.2f\n", rng[1], rng[2]))
if(rng[2] > 50) { expr <- log2(expr+1); cat("applied log2(x+1)\n") }
expr <- expr[complete.cases(expr),]

## --- limma DEG (batch-adjusted): Tumor vs Control ---
group <- relevel(group, ref=grep("ontrol", levels(group), value=TRUE)[1])
design <- if(!is.null(batch)) model.matrix(~0+group+batch) else model.matrix(~0+group)
gl <- make.names(levels(group)); colnames(design)[1:length(gl)] <- gl
tumor_lvl <- grep("umor", gl, value=TRUE)[1]; ctrl_lvl <- grep("ontrol", gl, value=TRUE)[1]
fit <- lmFit(expr, design)
cm <- makeContrasts(contrasts=paste0(tumor_lvl,"-",ctrl_lvl), levels=design)
fit2 <- eBayes(contrasts.fit(fit, cm), trend=TRUE, robust=TRUE)
tt <- topTable(fit2, number=Inf, sort.by="P")
nfdr <- sum(tt$adj.P.Val<0.05); nfdr_fc <- sum(tt$adj.P.Val<0.05 & abs(tt$logFC)>1)
nobs <- sum(tt$P.Value<0.05)
cat(sprintf("DEG FDR<0.05: %d | FDR<0.05 & |logFC|>1: %d | nominal p<0.05: %d (of %d)\n",
            nfdr, nfdr_fc, nobs, nrow(tt)))

## --- permutation signal test (shuffle group within batch if batch present) ---
nperm <- 500; permc <- integer(nperm)
for(i in seq_len(nperm)){
  if(!is.null(batch)){
    gp <- group
    for(b in levels(batch)){ idx <- which(batch==b); gp[idx] <- sample(group[idx]) }
  } else gp <- factor(sample(as.character(group)), levels=levels(group))
  d <- if(!is.null(batch)) model.matrix(~0+gp+batch) else model.matrix(~0+gp)
  gln <- make.names(levels(gp)); colnames(d)[1:length(gln)] <- gln
  f <- tryCatch(eBayes(contrasts.fit(lmFit(expr,d),
        makeContrasts(contrasts=paste0(grep("umor",gln,value=TRUE)[1],"-",grep("ontrol",gln,value=TRUE)[1]),levels=d)),
        trend=TRUE, robust=TRUE), error=function(e) NULL)
  permc[i] <- if(is.null(f)) NA else sum(f$p.value[,1]<0.05)
}
permc <- permc[!is.na(permc)]
emp_p <- (sum(permc>=nobs)+1)/(length(permc)+1)
cat(sprintf("PERMUTATION: observed p<0.05=%d | null mean=%.0f sd=%.0f max=%.0f | empirical p=%.4f\n",
            nobs, mean(permc), sd(permc), max(permc), emp_p))

fwrite(data.table(probe=rownames(tt), tt[,c("logFC","P.Value","adj.P.Val")]), file.path(outdir,"deg_probe_gate.csv"))
verdict <- if(nfdr>=300 && emp_p<0.01) "GO (strong real signal)" else if(nfdr>=50) "BORDERLINE" else "NO-GO"
cat(sprintf("\n=== SIGNAL GATE VERDICT: %s ===\n", verdict))
cat(sprintf("summary: FDR<0.05=%d, FDR&|FC|>1=%d, perm emp.p=%.4f\n", nfdr, nfdr_fc, emp_p))
writeLines(c(sprintf("verdict: %s", verdict),
             sprintf("DEG_FDR05: %d", nfdr),
             sprintf("DEG_FDR05_FC1: %d", nfdr_fc),
             sprintf("nominal_p05: %d", nobs),
             sprintf("perm_emp_p: %.4f", emp_p),
             sprintf("n_samples: %d", ncol(expr)),
             sprintf("group_counts: %s", paste(names(table(group)),table(group),collapse="; "))),
           file.path(outdir,"signal_gate_verdict.txt"))
cat("DONE\n")
