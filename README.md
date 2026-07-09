# Vestibular Schwannoma Composition-Aware Transcriptomics

This staging directory is prepared for the public GitHub repository:

`https://github.com/1210loku-lab/vestibular-schwannoma-composition-transcriptomics`

It accompanies the manuscript:

**Cross-cohort and single-cell decomposition of reproducible transcriptomic features in vestibular schwannoma**

This is a confirmatory, composition-aware reanalysis of public vestibular schwannoma transcriptomic datasets. It evaluates cross-platform reproducibility of bulk VS-versus-nerve expression features and maps reproducible signals onto single-cell-defined compartments using pseudobulk analysis, reference-based deconvolution and conservative composition-adjusted modelling.

## Public Datasets

- GSE39645: discovery bulk microarray cohort.
- GSE108524: external bulk microarray cohort.
- GSE230375: single-cell reference with VS and great-auricular-nerve samples.
- GSE216783: VS-only single-cell cohort used for reference-sensitivity analysis.

Raw public data should be obtained from NCBI GEO. This repository does not redistribute raw GEO downloads, single-cell objects, local working data, submission files, or author/internal handoff notes.

## Repository Contents

- `scripts/`: main analysis and submission-asset scripts.
- `source_data/`: figure source-data tables exported from the analysis outputs.
- `supplementary_tables/`: Supplementary Table S1, the compartment-aware cross-evidence map.
- `session_info/`: session information available for the single-cell localisation run.
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

## Reproducibility Notes

The analyses use public datasets with heterogeneous nerve controls and a limited number of single-cell nerve donors. Reference-based deconvolution values should be interpreted as relative computational estimates rather than measured tissue proportions. Candidate genes are provided as compartment-aware follow-up examples, not as diagnostic biomarkers or therapeutic targets.

The public-release scripts use `Sys.getenv("PROJ_ROOT", getwd())` for project-root resolution. Before running from a fresh clone, download the public GEO datasets, recreate the expected `data/raw/` and `data/processed/` inputs, then either run commands from the repository root or set `PROJ_ROOT=/path/to/clone`.

## Pre-Upload Leak Checks

Before creating the GitHub repository or Zenodo archive, run leak checks from this staging directory:

```bash
find . -type f -size +10M -print
git init
git config user.name "[Author/Lab]"
git config user.email "noreply@example.com"
git status --short
git remote -v
```

Also run a local path and secret-token scan before upload, including absolute user paths, local-file URL schemes, private-key headers, GitHub tokens, cloud access keys, Slack tokens, personal email addresses and cross-project keywords. Expected pre-upload status: no files above 10 MB, no secret/path hits, no remote until the user explicitly approves repository creation.

After the GitHub repository is created but before making it public, perform a fresh clone and repeat both file-tree and full-history checks. The 07 project release showed that a later commit can reintroduce an absolute path even after the initial Zenodo archive is clean, so the final check must cover current files, commit authors/committers and full history.

For Zenodo, cite the concept DOI in the manuscript and portal metadata once the repository release is archived. Version DOIs may change with later releases.

## License

Code is prepared under the MIT License. Public datasets remain governed by their original repositories and terms of use.
