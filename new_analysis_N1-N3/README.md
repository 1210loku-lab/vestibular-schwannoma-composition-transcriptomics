# New analysis code and source data

This directory contains the complete scripts and compact outputs for the three analyses added for the Genomics submission.

- `N1_permutation_null_concordance.R`: 10,000-permutation null for cross-cohort direction concordance and Pearson correlation.
- `N1_permutation_null_summary.csv`: observed and null summary statistics.
- `N2_three_way_ablation.R`: naive-versus-composition-adjusted threshold reclassification and compartment annotation.
- `N2_ablation_summary.csv`, `N2_classification_reproduce_check.csv`, `N2_compartment_by_class.csv`: reclassification and independent consistency checks.
- `N3_evidence_tier_resource.R`: construction and validation of the evidence-tiered 354-transcript resource.
- `N3_candidate_evidence_tier_resource.csv`, `N3_tier1_candidates.csv`, `N3_tier_summary.csv`: complete resource and tier summaries.

The scripts use project-relative source files documented in their headers. The full analysis archive and session-information records are available at https://doi.org/10.5281/zenodo.21361686. These additions are prepared for inclusion in the next version of that archive.
