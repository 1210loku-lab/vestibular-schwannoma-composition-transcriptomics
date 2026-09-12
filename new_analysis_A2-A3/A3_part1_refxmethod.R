#!/usr/bin/env Rscript
# A3 Part1: reference x method robustness matrix (isolate reference effect vs method effect)
options(stringsAsFactors=FALSE, warn=1)
suppressPackageStartupMessages({library(dtangle); library(nnls)})
ROOT<-normalizePath(Sys.getenv("PROJ_ROOT", getwd()), mustWork=TRUE)
OUT<-file.path(ROOT,"results","deconv_multimethod"); dir.create(OUT,showWarnings=FALSE)
CT<-c("Schwann/tumour","Myeloid","T cell","Fibroblast","Endothelial","Pericyte")  # canonical order

read_expr<-function(p){x<-read.csv(p,check.names=FALSE);g<-x[[1]];m<-as.matrix(x[,-1,drop=FALSE]);storage.mode(m)<-"double";rownames(m)<-g;m}
to_linear<-function(e){x<-as.matrix(e);storage.mode(x)<-"double";sh<-min(x,na.rm=TRUE);y<-2^(x-sh);y[!is.finite(y)|y<=0]<-.Machine$double.eps;dimnames(y)<-dimnames(x);y}
sel_markers<-function(sig_lin,n=150){ # sig_lin genes x celltypes (linear)
  mx<-sapply(seq_len(ncol(sig_lin)),function(j) apply(sig_lin[,-j,drop=FALSE],1,max))
  spec<-log2((sig_lin+.05)/(mx+.05))
  unique(unlist(lapply(seq_len(ncol(sig_lin)),function(j){e<-which(sig_lin[,j]>.05&is.finite(spec[,j]));head(rownames(sig_lin)[e[order(spec[e,j],decreasing=TRUE)]],n)})))}
run_nnls<-function(ref,bulk){com<-intersect(rownames(ref),rownames(bulk));ref<-as.matrix(ref[com,]);bul<-as.matrix(bulk[com,])
  k<-rowSums(is.finite(ref))==ncol(ref)&apply(bul,1,function(z)all(is.finite(z)))&rowSums(ref)>0;ref<-ref[k,];bul<-bul[k,]
  est<-t(vapply(seq_len(ncol(bul)),function(i){p<-pmax(coef(nnls(ref,bul[,i])),0);if(!sum(p)>0)p<-rep(1,length(p));p/sum(p)},numeric(ncol(ref))))
  rownames(est)<-colnames(bul);colnames(est)<-colnames(ref);est}
run_dtangle<-function(sig_lin,bulk_lin,mk_genes){
  com<-Reduce(intersect,list(rownames(bulk_lin),rownames(sig_lin),mk_genes))
  Y<-t(log2(bulk_lin[com,,drop=FALSE]+1)); ref<-t(log2(sig_lin[com,,drop=FALSE]+1))
  mk<-lapply(seq_len(nrow(ref)),function(i){s<-ref[i,]-apply(ref[-i,,drop=FALSE],2,max);order(s,decreasing=TRUE)[1:min(100,length(s))]})
  est<-dtangle(Y,reference=ref,markers=mk,data_type="microarray-gene")$estimates;colnames(est)<-rownames(ref);est}

# signatures -> linear, canonical column order
sig230<-read_expr(file.path(ROOT,"data","processed","GSE230375_celltype_meanexpr.csv")); colnames(sig230)<-CT; sig230<-expm1(pmax(sig230,0))
sig216<-read_expr(file.path(ROOT,"results","gse216783","celltype_meanexpr_gse216783.csv")); colnames(sig216)<-CT; sig216<-expm1(pmax(sig216,0))
refs<-list(GSE230375=sig230, GSE216783=sig216)
cohorts<-c("GSE39645","GSE108524")

long<-list()
for(co in cohorts){
  expr<-read_expr(file.path(ROOT,"data","processed",paste0(co,"_expr_gene.csv"))); bulk<-to_linear(expr)
  meta<-read.csv(file.path(ROOT,"data","processed",paste0(co,"_meta.csv"))); if(co=="GSE108524") meta$group<-ifelse(meta$group=="Control","Control","Tumor")
  grp<-meta$group[match(colnames(bulk),meta$sample)]
  for(rf in names(refs)){
    sig<-refs[[rf]]; mk<-sel_markers(sig)
    nnls_est<-run_nnls(sig[mk,,drop=FALSE],bulk)[,CT,drop=FALSE]
    dt_est<-run_dtangle(sig,bulk,mk)[,CT,drop=FALSE]
    for(mth in c("NNLS","dtangle")){
      est<-if(mth=="NNLS") nnls_est else dt_est
      for(ct in CT) long[[paste(co,rf,mth,ct)]]<-data.frame(cohort=co,reference=rf,method=mth,cell_type=ct,
        sample=rownames(est),group=grp[match(rownames(est),colnames(bulk))],fraction=est[,ct])
    }
  }
  # fold in existing MuSiC(GSE230375 ref)
  mf<-file.path(ROOT,"results","deconv",paste0("cell_fractions_",co,"_music.csv"))
  if(file.exists(mf)){x<-read.csv(mf,check.names=FALSE);cc<-intersect(colnames(x),c("Schwann/肿瘤细胞","巨噬/髓系","T细胞","成纤维","内皮","周细胞"))
    if(length(cc)==6){m<-as.matrix(x[,cc]);colnames(m)<-CT;rownames(m)<-x$sample
      for(ct in CT) long[[paste(co,"GSE230375","MuSiC",ct)]]<-data.frame(cohort=co,reference="GSE230375",method="MuSiC",cell_type=ct,sample=x$sample,group=x$group,fraction=m[,ct])}}
}
L<-do.call(rbind,long); write.csv(L,file.path(OUT,"reference_by_method_fractions_long.csv"),row.names=FALSE)

# stability: method effect (within reference, across methods) & reference effect (within method, across references)
sp<-function(a,b){if(sd(a,na.rm=TRUE)==0||sd(b,na.rm=TRUE)==0) return(NA); suppressWarnings(cor(a,b,method="spearman",use="complete.obs"))}
wide<-function(sub){ reshape(sub[,c("sample","fraction")],direction="wide",idvar="sample",timevar=NULL)}
rob<-list()
for(co in cohorts) for(ct in CT){
  d<-L[L$cohort==co&L$cell_type==ct,]
  cfg<-split(d$fraction,paste(d$reference,d$method)); # align by sample
  al<-lapply(split(d,paste(d$reference,d$method)),function(z)setNames(z$fraction,z$sample))
  samp<-Reduce(intersect,lapply(al,names)); al<-lapply(al,function(v)v[samp])
  # method effect within each ref
  me<-c(); for(rf in unique(d$reference)){k<-grep(paste0("^",rf," "),names(al),value=TRUE); if(length(k)>=2){cb<-combn(k,2); me<-c(me,apply(cb,2,function(p)sp(al[[p[1]]],al[[p[2]]])))}}
  # reference effect within each method
  re<-c(); for(mth in unique(d$method)){k<-grep(paste0(" ",mth,"$"),names(al),value=TRUE); if(length(k)>=2){cb<-combn(k,2); re<-c(re,apply(cb,2,function(p)sp(al[[p[1]]],al[[p[2]]])))}}
  rob[[paste(co,ct)]]<-data.frame(cohort=co,cell_type=ct,n_config=length(al),
    method_effect_min=ifelse(length(me),min(me,na.rm=TRUE),NA),method_effect_med=ifelse(length(me),median(me,na.rm=TRUE),NA),
    reference_effect_min=ifelse(length(re),min(re,na.rm=TRUE),NA),reference_effect_med=ifelse(length(re),median(re,na.rm=TRUE),NA))
}
R<-do.call(rbind,rob); write.csv(R,file.path(OUT,"reference_by_method_stability.csv"),row.names=FALSE)

# group contrasts per config with FDR
gc<-list()
for(co in cohorts) for(rf in unique(L$reference[L$cohort==co])) for(mth in unique(L$method[L$cohort==co&L$reference==rf])){
  for(ct in CT){d<-L[L$cohort==co&L$reference==rf&L$method==mth&L$cell_type==ct,]
    tu<-d$fraction[d$group=="Tumor"];co2<-d$fraction[d$group=="Control"]
    p<-tryCatch(wilcox.test(tu,co2,exact=FALSE)$p.value,error=function(e)NA)
    gc[[paste(co,rf,mth,ct)]]<-data.frame(cohort=co,reference=rf,method=mth,cell_type=ct,n_tumor=length(tu),n_ctrl=length(co2),median_shift=median(tu)-median(co2),p=p)}}
G<-do.call(rbind,gc); G$FDR<-ave(G$p,paste(G$cohort,G$reference,G$method),FUN=function(x)p.adjust(x,"BH"))
write.csv(G,file.path(OUT,"reference_by_method_group_contrasts.csv"),row.names=FALSE)

cat("=== stability (Spearman) ===\n"); print(R,row.names=FALSE,digits=2)
cat("\n=== group direction agreement across ALL configs per cohort x cell_type ===\n")
for(co in cohorts) for(ct in CT){s<-sign(G$median_shift[G$cohort==co&G$cell_type==ct]);a<-length(unique(s[!is.na(s)]))==1
  cat(sprintf("%-10s %-15s configs=%d all_same_dir=%s anySig(FDR<.05)=%s\n",co,ct,sum(G$cohort==co&G$cell_type==ct),a,any(G$FDR[G$cohort==co&G$cell_type==ct]<0.05,na.rm=TRUE)))}
cat("\nsaved:",OUT,"\n")
