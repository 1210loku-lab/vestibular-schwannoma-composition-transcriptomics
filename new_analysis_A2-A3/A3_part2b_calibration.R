#!/usr/bin/env Rscript
# A3 Part2b: re-run known-mixture with per-mixture detail + MAE/CCC/calibration slope + donor CI.
options(stringsAsFactors=FALSE, warn=1); set.seed(1)
suppressPackageStartupMessages({library(data.table); library(dtangle); library(nnls)})
ROOT<-normalizePath(Sys.getenv("PROJ_ROOT", getwd()), mustWork=TRUE)
OUT<-file.path(ROOT,"results","deconv_multimethod")
CT<-c("Schwann/tumour","Myeloid","T cell","Fibroblast","Endothelial","Pericyte")
CNT<-file.path(ROOT,"data","raw","GSE216783_sch.all.counts.tsv.gz")
meta<-fread(file.path(ROOT,"results","gse216783","cell_meta.csv")); meta<-meta[cell_type %in% CT]
cap<-200; sub<-meta[, .SD[sample(.N,min(.N,cap))], by=.(sample,cell_type)]
tab<-dcast(sub[, .N, by=.(sample,cell_type)], sample~cell_type, value.var="N", fill=0)
usable<-tab$sample[apply(as.matrix(tab[,-1,with=FALSE]),1,function(r) all(r>=15))]
keepbc<-sub$barcode; hdr<-strsplit(readLines(CNT,n=1),"\t")[[1]]; sel<-which(hdr %in% keepbc)
dt<-fread(CNT, select=c(1, sel+1), header=FALSE, skip=1); genes<-dt[[1]]; M<-as.matrix(dt[,-1,with=FALSE]); rownames(M)<-genes; colnames(M)<-hdr[sel]
M<-M[, intersect(keepbc,colnames(M)),drop=FALSE]; M<-M[rowSums(M)>0,]
bc2ct<-setNames(sub$cell_type,sub$barcode); bc2do<-setNames(sub$sample,sub$barcode)
cellct<-bc2ct[colnames(M)]; celldo<-bc2do[colnames(M)]
cpm<-function(m){s<-colSums(m);s[s==0]<-1;sweep(m,2,s,"/")*1e6}
sel_markers<-function(sig,n=150){mx<-sapply(seq_len(ncol(sig)),function(j)apply(sig[,-j,drop=FALSE],1,max));spec<-log2((sig+.05)/(mx+.05));unique(unlist(lapply(seq_len(ncol(sig)),function(j){e<-which(sig[,j]>.05&is.finite(spec[,j]));head(rownames(sig)[e[order(spec[e,j],decreasing=TRUE)]],n)})))}
nnls_one<-function(ref,y){com<-intersect(rownames(ref),names(y));r<-ref[com,,drop=FALSE];yy<-y[com];k<-rowSums(r)>0&is.finite(yy);p<-pmax(coef(nnls(as.matrix(r[k,]),yy[k])),0);if(!sum(p)>0)p<-rep(1,length(p));setNames(p/sum(p),colnames(ref))}
dtangle_bulk<-function(sig,mk,Ymat){com<-Reduce(intersect,list(rownames(Ymat),rownames(sig),mk));Y<-t(log2(Ymat[com,,drop=FALSE]+1));ref<-t(log2(sig[com,,drop=FALSE]+1));mkl<-lapply(seq_len(nrow(ref)),function(i){s<-ref[i,]-apply(ref[-i,,drop=FALSE],2,max);order(s,decreasing=TRUE)[1:min(100,length(s))]});e<-dtangle(Y,reference=ref,markers=mkl,data_type="rna-seq")$estimates;colnames(e)<-rownames(ref);e}
sig_from<-function(cols){cc<-colnames(M)[cols];ct<-cellct[cc];sapply(CT,function(t){idx<-cc[ct==t];if(length(idx)==0)return(rep(0,nrow(M)));rowMeans(cpm(M[,idx,drop=FALSE]))})}
Mmix<-100; ncell<-500; detail<-list()
for(hd in usable){
  tr<-which(celldo!=hd); te<-which(celldo==hd); sig<-sig_from(tr); sig<-sig[rowSums(sig)>0,]; mk<-sel_markers(sig)
  tb<-colnames(M)[te]; tct<-cellct[tb]
  props<-matrix(rgamma(Mmix*6,1),Mmix,6); props<-props/rowSums(props); colnames(props)<-CT
  Ymat<-matrix(0,nrow(M),Mmix,dimnames=list(rownames(M),NULL)); truep<-matrix(0,Mmix,6,dimnames=list(NULL,CT))
  for(i in 1:Mmix){n_i<-round(props[i,]*ncell);n_i[n_i<1&props[i,]*ncell>0]<-1
    pk<-unlist(lapply(CT,function(t){pool<-tb[tct==t];if(length(pool)==0||n_i[t]==0)return(character(0));sample(pool,n_i[t],replace=length(pool)<n_i[t])}))
    truep[i,]<-table(factor(cellct[pk],levels=CT))/length(pk); Ymat[,i]<-rowSums(M[,pk,drop=FALSE])}
  Ymat<-cpm(Ymat)
  nn<-t(sapply(1:Mmix,function(i)nnls_one(sig[mk,,drop=FALSE],Ymat[,i])))[,CT,drop=FALSE]
  dtn<-tryCatch(dtangle_bulk(sig,mk,Ymat)[,CT,drop=FALSE],error=function(e)matrix(NA,Mmix,6,dimnames=list(NULL,CT)))
  for(mth in c("NNLS","dtangle")){est<-if(mth=="NNLS")nn else dtn
    for(t in CT) detail[[paste(hd,mth,t)]]<-data.frame(donor=hd,method=mth,cell_type=t,mix=1:Mmix,true=truep[,t],est=est[,t])}
  cat("done",hd,"\n")
}
D<-do.call(rbind,detail); fwrite(D,file.path(OUT,"known_mixture_permixture.csv"))
ccc<-function(x,y){mx<-mean(x);my<-mean(y);sx<-var(x);sy<-var(y);sxy<-cov(x,y);2*sxy/(sx+sy+(mx-my)^2)}
# per method x cell_type (pooled over donors/mixtures) + donor-level CI on pearson
agg<-list()
for(mth in c("NNLS","dtangle")) for(t in CT){d<-D[D$method==mth&D$cell_type==t,]
  slope<-coef(lm(est~true,d))[2]; per_donor_r<-sapply(split(d,d$donor),function(z)suppressWarnings(cor(z$true,z$est)))
  agg[[paste(mth,t)]]<-data.frame(method=mth,cell_type=t,
    pearson=cor(d$true,d$est),CCC=ccc(d$true,d$est),MAE=mean(abs(d$est-d$true)),RMSE=sqrt(mean((d$est-d$true)^2)),
    bias=mean(d$est-d$true),calib_slope=slope,
    donor_r_min=min(per_donor_r,na.rm=TRUE),donor_r_lo=quantile(per_donor_r,.025,na.rm=TRUE),donor_r_hi=quantile(per_donor_r,.975,na.rm=TRUE))}
A<-do.call(rbind,agg); fwrite(A,file.path(OUT,"known_mixture_calibration.csv"))
cat("\n=== calibration (pooled) ===\n"); print(A,row.names=FALSE,digits=2)
# rare-proportion subset (true<0.10)
rare<-D[D$true<0.10,]; ra<-list()
for(mth in c("NNLS","dtangle")) for(t in CT){d<-rare[rare$method==mth&rare$cell_type==t,]; if(nrow(d)<5)next
  ra[[paste(mth,t)]]<-data.frame(method=mth,cell_type=t,n=nrow(d),MAE_rare=mean(abs(d$est-d$true)),bias_rare=mean(d$est-d$true))}
RA<-do.call(rbind,ra); fwrite(RA,file.path(OUT,"known_mixture_rare_lowprop.csv"))
cat("\n=== rare (true<0.10) MAE/bias ===\n"); print(RA,row.names=FALSE,digits=3)
cat("saved per-mixture + calibration + rare tables\n")
