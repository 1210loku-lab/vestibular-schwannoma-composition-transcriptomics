import warnings; warnings.filterwarnings("ignore")
import scanpy as sc, numpy as np, pandas as pd, anndata
sc.settings.verbosity=0
np.random.seed(0)

files={'AN014':'gse/GSM7973448_AN014_filtered_feature_bc_matrix.h5',
       'AN017':'gse/GSM7973449_AN017_filtered_feature_bc_matrix.h5',
       'AN018':'gse/GSM7973450_AN018_filtered_feature_bc_matrix.h5'}
ads=[]
for s,fn in files.items():
    a=sc.read_10x_h5(fn); a.var_names_make_unique(); a.obs['donor']=s; ads.append(a)
adata=anndata.concat(ads,label='sample',keys=list(files),index_unique='-',join='outer')
adata.var['mt']=adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata,qc_vars=['mt'],inplace=True,percent_top=None)
# strict QC matching original paper: 200<nFeature<6000, mt<10
adata=adata[(adata.obs.n_genes_by_counts>200)&(adata.obs.n_genes_by_counts<6000)&(adata.obs.pct_counts_mt<10)].copy()
sc.pp.filter_genes(adata,min_cells=3)
print("post-QC:",adata.shape, "| per-donor:",dict(adata.obs.donor.value_counts()))
adata.layers['counts']=adata.X.copy()
sc.pp.normalize_total(adata,target_sum=1e4); sc.pp.log1p(adata); adata.raw=adata
sc.pp.highly_variable_genes(adata,n_top_genes=2000,batch_key='donor')
ad=adata[:,adata.var.highly_variable].copy()
sc.pp.scale(ad,max_value=10); sc.tl.pca(ad,n_comps=30)
# harmony integration across donors
sc.external.pp.harmony_integrate(ad,'donor',max_iter_harmony=20)
sc.pp.neighbors(ad,n_neighbors=15,n_pcs=30,use_rep='X_pca_harmony')
sc.tl.leiden(ad,resolution=1.0,flavor='igraph',n_iterations=2,directed=False)
sc.tl.umap(ad)
adata.obs['leiden']=ad.obs['leiden'].values; adata.obsm['X_umap']=ad.obsm['X_umap']
print("clusters:",adata.obs['leiden'].nunique())

markers={'Schwann/tumour':['SOX10','S100B','PLP1','MPZ','NGFR','SOX2'],
         'Myeloid':['PTPRC','CD68','LYZ','AIF1','C1QB','ITGAM'],
         'T cell':['CD3E','CD3D','TRAC','CD2'],
         'Endothelial':['PECAM1','CLDN5','VWF','FLT1','CDH5'],
         'Fibroblast':['DCN','COL1A1','COL1A2','LUM','COL3A1'],
         'Pericyte':['RGS5','ACTA2','PDGFRB','NOTCH3','MCAM']}
for lin,gs in markers.items():
    g=[x for x in gs if x in adata.raw.var_names]
    sc.tl.score_genes(adata,g,score_name='sc_'+lin,use_raw=True)
scorecols=['sc_'+l for l in markers]
cl_ann={}
for cl in sorted(adata.obs['leiden'].unique(),key=int):
    m=(adata.obs['leiden']==cl).values
    means={l:adata.obs.loc[m,'sc_'+l].mean() for l in markers}
    cl_ann[cl]=max(means,key=means.get)
adata.obs['celltype']=adata.obs['leiden'].map(cl_ann).astype(str)
print("\ncluster->celltype:",cl_ann)
print("celltype sizes:\n",adata.obs['celltype'].value_counts().to_string())
fibro_present = 'Fibroblast' in set(cl_ann.values())
print("\n*** Fibroblast cluster recovered? ->", fibro_present, "***")
# show fibroblast marker top cluster regardless
print("Mean DCN/COL1A1 by cluster (top5):")
for g in ['DCN','COL1A1']:
    if g in adata.raw.var_names:
        x=np.asarray(adata.raw[:,g].X.todense()).ravel()
        s=pd.Series(x,index=adata.obs_names).groupby(adata.obs['leiden']).mean().sort_values(ascending=False)
        print(f"  {g}: "+", ".join(f"cl{k}={v:.2f}({cl_ann[k]})" for k,v in s.head(5).items()))

# candidate localization, FULL denominator, donor-aware, specificity margin
pm=pd.read_csv('candidate_priority_matrix.csv')[['gene','pb_dominant_celltype']].dropna(subset=['gene'])
cts=list(markers)
rows=[]
for _,r in pm.iterrows():
    g=r['gene']; assigned=r['pb_dominant_celltype']
    if g not in adata.raw.var_names:
        rows.append([g,assigned,'GENE_ABSENT','',0,0,0]); continue
    x=np.asarray(adata.raw[:,g].X.todense()).ravel()
    mean_by={ct:float(x[(adata.obs['celltype']==ct).values].mean()) if (adata.obs['celltype']==ct).sum()>0 else 0 for ct in cts}
    ordered=sorted(mean_by,key=mean_by.get,reverse=True)
    dom=ordered[0]; second=ordered[1]
    margin=mean_by[dom]-mean_by[second]
    # donor-aware detection
    ndon=sum(1 for d in files if (x[(adata.obs['donor']==d).values]>0).mean()>0.01)
    conc='YES' if dom==assigned else ('PARTIAL' if assigned in ordered[:2] else 'NO')
    rows.append([g,assigned,dom,conc,round(100*float((x>0).mean()),2),round(margin,3),ndon])
res=pd.DataFrame(rows,columns=['gene','manuscript_compartment','scRNA_dominant','concordance','pct_nonzero','spec_margin_log','n_donors_detected'])
res.to_csv('A2_redo_localization.csv',index=False)
nonabs=res[res.scRNA_dominant!='GENE_ABSENT']
def rate(d):
    n=len(d); return n,(d.concordance=='YES').sum(),(d.concordance.isin(['YES','PARTIAL'])).sum()
n,y,t=rate(nonabs)
print(f"\n=== FULL-denominator concordance (n={n}) ===  exact {y}={100*y/n:.1f}%  top2 {t}={100*t/n:.1f}%")
# confident subset: dominant detected in >=2 donors
conf=nonabs[nonabs.n_donors_detected>=2]
nc,yc,tc=rate(conf)
print(f"Confident subset (dominant gene in >=2 donors, n={nc}): exact {100*yc/nc:.1f}%  top2 {100*tc/nc:.1f}%")
print("\nBy manuscript compartment (full):")
for ct,gp in nonabs.groupby('manuscript_compartment'):
    nn,yy,tt=rate(gp); print(f"  {ct:16s} n={nn:3d}  exact {100*yy/nn:4.0f}%  top2 {100*tt/nn:4.0f}%")
print("\nNamed 4:")
print(nonabs[nonabs.gene.isin(['L1CAM','NRG1','FCGBP','PRRX1'])][['gene','manuscript_compartment','scRNA_dominant','concordance','pct_nonzero','n_donors_detected']].to_string(index=False))

import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sc.pl.umap(adata,color=['celltype'],show=False); plt.tight_layout(); plt.savefig('A2_redo_umap.png',dpi=150,bbox_inches='tight'); plt.close()
g4=[g for g in ['L1CAM','NRG1','FCGBP','PRRX1'] if g in adata.raw.var_names]
sc.pl.dotplot(adata,g4,groupby='celltype',show=False); plt.savefig('A2_redo_dotplot.png',dpi=150,bbox_inches='tight'); plt.close()
print("\nsaved A2_redo_localization.csv, A2_redo_umap.png, A2_redo_dotplot.png")
