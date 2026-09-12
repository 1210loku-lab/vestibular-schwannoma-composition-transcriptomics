#!/usr/bin/env python3
"""A2 independent-cohort compartment concordance: point estimate, bootstrap 95% CI,
label-permutation chance baseline and empirical p, and confusion matrix.
Input: A2_redo_localization.csv (per-transcript manuscript vs GSE250061 dominant compartment).
Reproducible: fixed seed, 10,000 resamples/permutations."""
import pandas as pd, numpy as np
np.random.seed(0)
r = pd.read_csv('A2_redo_localization.csv')
na = r[r.scRNA_dominant != 'GENE_ABSENT'].copy()
n = len(na)
exact = (na.concordance == 'YES').mean()
top2  = (na.concordance.isin(['YES','PARTIAL'])).mean()
man, dom = na.manuscript_compartment.values, na.scRNA_dominant.values
# label-permutation chance baseline + empirical p (exact match)
perm = np.array([np.mean(np.random.permutation(man) == dom) for _ in range(10000)])
chance = perm.mean(); pval = (np.sum(perm >= exact) + 1) / (len(perm) + 1)
# bootstrap 95% CI on exact
idx = np.arange(n)
bs = np.array([(na.concordance.values[np.random.choice(idx, n, replace=True)] == 'YES').mean() for _ in range(10000)])
lo, hi = np.percentile(bs, [2.5, 97.5])
with open('A2_concordance_stats_output.txt', 'w') as f:
    f.write(f"n evaluable transcripts: {n}\n")
    f.write(f"exact concordance: {exact*100:.1f}% (95% bootstrap CI {lo*100:.1f}-{hi*100:.1f}%)\n")
    f.write(f"top-2 concordance: {top2*100:.1f}%\n")
    f.write(f"label-permutation chance (exact): {chance*100:.1f}%; empirical p = {pval:.4g} (10,000 permutations)\n")
    by = na.groupby('manuscript_compartment').apply(
        lambda d: pd.Series({'n': len(d), 'exact_%': round(100*(d.concordance=='YES').mean(),1),
                             'top2_%': round(100*d.concordance.isin(['YES','PARTIAL']).mean(),1)}))
    f.write("\nby manuscript compartment:\n" + by.to_string() + "\n")
    cm = pd.crosstab(na.manuscript_compartment, na.scRNA_dominant)
    cm.to_csv('A2_confusion_matrix.csv')
    f.write("\nconfusion matrix saved to A2_confusion_matrix.csv\n")
print(open('A2_concordance_stats_output.txt').read())
