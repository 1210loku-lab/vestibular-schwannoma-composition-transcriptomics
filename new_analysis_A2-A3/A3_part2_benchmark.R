#!/usr/bin/env Rscript
# A3 Part2: donor-held-out known-mixture (pseudobulk) benchmark -> tests algorithmic
# recovery against known synthetic-mixture proportions.
options(stringsAsFactors=FALSE, warn=1); set.seed(1)
suppressPackageStartupMessages({library(data.table); library(dtangle); library(nnls)})
ROOT<-normalizePath(Sys.getenv("PROJ_ROOT", getwd()), mustWork=TRUE)
OUT<-file.path(ROOT,"results","deconv_multimethod")
CT<-c("Schwann/tumour","Myeloid","T cell","Fibroblast","Endothelial","Pericyte")
CNT<-file.path(ROOT,"data","raw","GSE216783_sch.all.counts.tsv.gz")

meta<-fread(file.path(ROOT,"results","gse216783","cell_meta.csv"))
meta<-meta[cell_type %in% CT]
# subsample per donor x celltype to bound memory
cap<-200
sub<-meta[, .SD[sample(.N, min(.N,cap))], by=.(sample,cell_type)]
# donors with all 6 types >=15 cells (usable as held-out ground-truth source)
tab<-dcast(sub[, .N, by=.(sample,cell_type)], sample~cell_type, value.var="N", fill=0)
usable<-tab$sample[apply(as.matrix(tab[,-1,with=FALSE]),1,function(r) all(r>=15))]
cat("usable held-out donors (all 6 types >=15 cells):",length(usable),"/",nrow(tab),"->",paste(usable,collapse=","),"\n")
keepbc<-sub$barcode

# read header to map barcode->col index
hdr<-strsplit(readLines(CNT,n=1),"\t")[[1]]           # length = n_cells (no gene col name)
# fread: first field per row is gene; columns after are cells in hdr order
sel_idx<-which(hdr %in% keepbc)
cat("reading",length(sel_idx),"cells x genes ...\n")
dt<-fread(CNT, select=c(1, sel_idx+1), header=FALSE, skip=1)   # col1=gene, +1 offset for gene col
genes<-dt[[1]]; M<-as.matrix(dt[,-1,with=FALSE]); rownames(M)<-genes
colnames(M)<-hdr[sel_idx]
M<-M[, intersect(keepbc, colnames(M)), drop=FALSE]
bc2ct<-setNames(sub$cell_type, sub$barcode); bc2do<-setNames(sub$sample, sub$barcode)
cellct<-bc2ct[colnames(M)]; celldo<-bc2do[colnames(M)]
cat("counts matrix:",nrow(M),"genes x",ncol(M),"cells\n")

cpm<-function(m){ s<-colSums(m); s[s==0]<-1; sweep(m,2,s,"/")*1e6 }
sig_from<-function(cols){ cc<-colnames(M)[cols]; ct<-cellct[cc]
  sapply(CT,function(t){ idx<-cc[ct==t]; if(length(idx)==0) return(rep(0,nrow(M))); rowMeans(cpm(M[,idx,drop=FALSE])) }) } # genes x CT (linear cpm)
sel_markers<-function(sig,n=150){ mx<-sapply(seq_len(ncol(sig)),function(j) apply(sig[,-j,drop=FALSE],1,max)); spec<-log2((sig+.05)/(mx+.05))
  unique(unlist(lapply(seq_len(ncol(sig)),function(j){e<-which(sig[,j]>.05&is.finite(spec[,j]));head(rownames(sig)[e[order(spec[e,j],decreasing=TRUE)]],n)}))) }
nnls_one<-function(ref,y){com<-intersect(rownames(ref),names(y));r<-ref[com,,drop=FALSE];yy<-y[com];k<-rowSums(r)>0&is.finite(yy);p<-pmax(coef(nnls(as.matrix(r[k,]),yy[k])),0);if(!sum(p)>0)p<-rep(1,length(p));setNames(p/sum(p),colnames(ref))}
dtangle_bulk<-function(sig,mk,Ymat){com<-Reduce(intersect,list(rownames(Ymat),rownames(sig),mk))
  Y<-t(log2(Ymat[com,,drop=FALSE]+1));ref<-t(log2(sig[com,,drop=FALSE]+1))
  mkl<-lapply(seq_len(nrow(ref)),function(i){s<-ref[i,]-apply(ref[-i,,drop=FALSE],2,max);order(s,decreasing=TRUE)[1:min(100,length(s))]})
  e<-dtangle(Y,reference=ref,markers=mkl,data_type="rna-seq")$estimates;colnames(e)<-rownames(ref);e}

M<-M[rowSums(M)>0,]
res<-list(); Mmix<-100; ncell<-500
for(hd in usable){
  train<-which(celldo!=hd); test<-which(celldo==hd)
  sig<-sig_from(train); sig<-sig[rowSums(sig)>0,]; mk<-sel_markers(sig)
  testbc<-colnames(M)[test]; testct<-cellct[testbc]
  # build Mmix pseudobulks
  props<-matrix(rgamma(Mmix*6,1),Mmix,6); props<-props/rowSums(props); colnames(props)<-CT
  Ymat<-matrix(0,nrow(M),Mmix,dimnames=list(rownames(M),paste0("mix",1:Mmix)))
  truep<-matrix(0,Mmix,6,dimnames=list(NULL,CT))
  for(i in 1:Mmix){ n_i<-round(props[i,]*ncell); n_i[n_i<1&props[i,]*ncell>0]<-1
    picked<-unlist(lapply(CT,function(t){pool<-testbc[testct==t]; if(length(pool)==0||n_i[t]==0) return(character(0)); sample(pool,n_i[t],replace=length(pool)<n_i[t])}))
    truep[i,]<-table(factor(cellct[picked],levels=CT))/length(picked)
    Ymat[,i]<-rowSums(M[,picked,drop=FALSE]) }
  Ymat<-cpm(Ymat)
  # NNLS
  nn<-t(sapply(1:Mmix,function(i) nnls_one(sig[mk,,drop=FALSE], Ymat[,i])))[,CT,drop=FALSE]
  # dtangle
  dt<-tryCatch(dtangle_bulk(sig,mk,Ymat)[,CT,drop=FALSE],error=function(e)matrix(NA,Mmix,6,dimnames=list(NULL,CT)))
  for(mth in c("NNLS","dtangle")){ est<-if(mth=="NNLS") nn else dt
    for(t in CT){ tv<-truep[,t]; ev<-est[,t]
      res[[paste(hd,mth,t)]]<-data.frame(donor=hd,method=mth,cell_type=t,
        pearson=suppressWarnings(cor(tv,ev)),spearman=suppressWarnings(cor(tv,ev,method="spearman")),
        rmse=sqrt(mean((tv-ev)^2)),bias=mean(ev-tv)) } }
  cat("done held-out",hd,"\n")
}
B<-do.call(rbind,res)
agg<-aggregate(cbind(pearson,spearman,rmse,bias)~method+cell_type,B,function(x)mean(x,na.rm=TRUE))
write.csv(B,file.path(OUT,"known_mixture_benchmark_perdonor.csv"),row.names=FALSE)
write.csv(agg,file.path(OUT,"known_mixture_benchmark_summary.csv"),row.names=FALSE)
cat("\n=== known-mixture recovery (mean over",length(usable),"held-out donors) ===\n"); print(agg,row.names=FALSE,digits=2)
writeLines(capture.output(sessionInfo()), file.path(OUT,"A3_sessionInfo.txt"))
cat("\nsaved:",OUT,"\n")
