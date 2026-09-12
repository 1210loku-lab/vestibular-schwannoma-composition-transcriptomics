# Vestibular Schwannoma Composition-Aware Transcriptomics

This repository accompanies the public GitHub project:

`https://github.com/1210loku-lab/vestibular-schwannoma-composition-transcriptomics`

It accompanies the manuscript:

**Cell-composition sensitivity of bulk transcriptomic inference in vestibular schwannoma: cross-cohort, multi-reference analysis with external single-cell localization**

This is a confirmatory, composition-aware reanalysis of public vestibular schwannoma transcriptomic datasets. It evaluates cross-platform reproducibility of bulk VS-versus-nerve expression features and maps reproducible signals onto single-cell-defined compartments using pseudobulk analysis, reference-based deconvolution and conservative composition-adjusted modelling.

## Public Datasets

- GSE39645: discovery bulk microarray cohort.
- GSE108524: external bulk microarray cohort.
- GSE230375: single-cell reference with VS and great-auricular-nerve samples.
- GSE216783: VS-only single-cell cohort used for reference-sensitivity and donor-held-out synthetic-mixture analyses.
- GSE250061: external tumour-only single-cell cohort used for descriptive broad-lineage localization.

Raw public data should be obtained from NCBI GEO. This repository does not redistribute raw GEO downloads, single-cell objects, local working data, submission files, or author/internal handoff notes.

## Repository Contents

- `scripts/`: main analysis and submission-asset scripts.
- `source_data/`: figure source-data tables exported from the analysis outputs.
- `supplementary_tables/`: Supplementary Table S1, the compartment-aware cross-evidence map.
- `session_info/`: session information available for the single-cell localisation run.
- `new_analysis_N1-N3/`: added in v1.1.0 — code and compact outputs for the three analyses supporting the revised submission (permutation null for cross-cohort concordance, naive-versus-composition-adjusted reclassification, and the evidence-tiered candidate resource).
- `new_analysis_A2-A3/`: added in v1.2.0 — code and compact outputs for external GSE250061 localization, the multi-reference and multi-estimator robustness matrix, the donor-held-out synthetic-mixture benchmark, calibration and rare-proportion summaries, Figure 5 source data, and session information.
- `RELEASE_MANIFEST.csv`: file inventory for this staging release.
- `.gitignore`: public-release guardrails modelled after the 07 project workflow.

## Main Analysis Scripts

- `scripts/01_load_signal_gate.R`: dataset loading and signal-quality gate.
- `scripts/02_annotate_deg_enrich.R`: differential expression and enrichment analysis.
- `scripts/03_immune_deconv.R`: immune/stromal score analysis.
- `scripts/04_wgcna_hub.R`: WGCNA modules.
- `scripts/05_hub.R`: hub-gene filtering.
- `scripts/07_validate_GSE108524.R`: cross-platform validation.
- `scripts/08_scrna_localize.R`: single-cell localisation.
- `scripts/09_nf2_vs_sporadic.R`: NF2/sporadic label audit retained for provenance.
- `scripts/12_pseudobulk_DE.R`: sample-level pseudobulk analysis.
- `scripts/13_bulk_deconv.R`: MuSiC and NNLS deconvolution.
- `scripts/14_gse216783_second_ref.R`: GSE216783 reference-sensitivity analysis.
- `scripts/15_candidate_matrix.R`: compartment-aware cross-evidence map.
- `scripts/18_figures_v2.R`: final submission figures.
- `scripts/19_module_preservation.R`: WGCNA module preservation.
- `scripts/20_prepare_submission_assets.R`: source-data export.
- `scripts/21_build_embedded_docx.py`: DOCX generation helper used locally; included for provenance.
- `new_analysis_N1-N3/N1_permutation_null_concordance.R`: 10,000-permutation null for cross-cohort direction concordance and correlation.
- `new_analysis_N1-N3/N2_three_way_ablation.R`: naive-versus-composition-adjusted threshold reclassification and compartment annotation.
- `new_analysis_N1-N3/N3_evidence_tier_resource.R`: construction and validation of the evidence-tiered candidate resource.

## Reproducibility Notes

The analyses use public datasets with heterogeneous nerve controls and a limited number of single-cell nerve donors. Reference-based deconvolution values should be interpreted as relative computational estimates rather than measured tissue proportions. Candidate genes are provided as compartment-aware follow-up examples, not as diagnostic biomarkers or therapeutic targets.

The public-release scripts use `Sys.getenv("PROJ_ROOT", getwd())` for project-root resolution. Before running from a fresh clone, download the public GEO datasets, recreate the expected `data/raw/` and `data/processed/` inputs, then either run commands from the repository root or set `PROJ_ROOT=/path/to/clone`.

## Archiving

Versioned releases are archived through the GitHub-Zenodo integration under Zenodo concept DOI [10.5281/zenodo.21361686](https://doi.org/10.5281/zenodo.21361686). Cite the version-specific DOI when referring to a particular release.

## License

Code is prepared under the MIT License. Public datasets remain governed by their original repositories and terms of use.
