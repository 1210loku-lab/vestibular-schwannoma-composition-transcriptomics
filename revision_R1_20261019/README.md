# Vestibular schwannoma R1 analysis code and processed source data

**Release version:** v1.3.0 (JMN R1 revision code and processed source data). This version is identified by the corresponding GitHub release tag. Zenodo assigns a version-specific DOI when it automatically archives the release; use that version DOI when citing this snapshot.

## Contents

- `scripts/`: 14 R scripts adapted from the R1 bulk-model, donor-reference, simulation, and figure-rendering workflows. Main-figure renderers write `Figure_1`, `Figure_2`, and `Figure_3`; the supplementary renderer writes S15–S18.
- `inputs/`: 22 frozen processed inputs and additional main-figure source tables listed in `INPUT_REQUIREMENTS.csv`, with their original `data/processed/` and `results/` relative layout retained beneath this directory. Raw GEO downloads and a Seurat object are not included.
- `outputs/`: frozen analysis tables, generated simulation prerequisites, fit/metric tables, validator snapshots, ten frozen Fig. 1–3 plotted-data tables, and package-relative integrity manifests. The bulk and donor result tables are retained alongside simulation inputs and outputs.
- `metadata/INPUT_SHA256.csv`: package-relative source names and packaged hashes for all 22 frozen input tables.
- `metadata/SOURCE_MAP.csv`: reader-facing source roles, package-relative file locations, hashes, and high-level transformation summaries; no source-checkout directory paths are included.
- `MANIFEST.csv`: byte count and SHA-256 for each package file except the manifest itself.
- `RUNTIME_VERSIONS.csv` and `DEPENDENCIES.md`: recorded runtime/package evidence and the script-level dependency summary.
- `LICENSE`: the existing project MIT license, applying to code. Upstream datasets and derived data remain subject to their source terms.

## Data sources and scope

The bundled processed inputs cover GSE39645, GSE108524, GSE230375, and GSE216783, together with project-derived label, deconvolution, robustness, validation, and WGCNA summary tables from those public-data workflows. Nine additional Fig. 1–3 source tables are included in their source-relative input paths; Figure 1 also reuses the already-bundled `deg_composition_adjusted.csv` table. Ten frozen panel plotting-data tables are included for comparison with renderer output. The donor-reference scripts read processed pseudobulk counts and cell metadata; they do not require a Seurat RDS. No Fig. 4 or provisional UMAP source is included.

The included outputs are frozen snapshots from the local R1 analyses. The controlled simulation records 20 replicates in each of the overlap and moderate/separation scenarios; the strong extension records 20 strong-scenario replicates. In the bulk-model results, the separately labeled tissue-only MuSiC PC1–PC5 sensitivity has group VIF 10.353, above the prespecified limit of 10, and fails the original numeric gate. The R1 source package’s historical strong-validation snapshot has 12 of 13 checks passing, with its original design-gate reporting check failure retained; that check records seven all-route gate failures and estimated-MuSiC VIF values of 5.374/8.322/13.070 (minimum/median/maximum). These are frozen historical results, not fresh validation results or a claim that the portable pipeline passes.

## Reproduction layout and order

Run from a separate writable copy of this package. Keep `VS_INPUT_ROOT` pointed at the package’s `inputs/` directory and `VS_OUTPUT_ROOT` at the writable package-copy root. Do not use the review archive itself as an execution target, because the scripts write derived files and manifests under `outputs/` and figures under `figures/`.

For example, after copying and entering the work copy:

```sh
package_root="$(pwd)"
export VS_INPUT_ROOT="$package_root/inputs"
export VS_OUTPUT_ROOT="$package_root"
Rscript scripts/bulk_model_effect_agreement.R
Rscript scripts/donor_reference_analysis.R
Rscript scripts/prepare_controlled_simulation.R
Rscript scripts/probe_all_gene_simulation.R
Rscript scripts/run_controlled_simulation.R
Rscript scripts/validate_controlled_simulation.R
Rscript scripts/prepare_strong_simulation.R
Rscript scripts/probe_all_gene_strong_simulation.R
Rscript scripts/run_strong_simulation.R
Rscript scripts/validate_strong_simulation.R
Rscript scripts/render_publication_figures.R
```

The donor-reference script sources `scripts/legacy_donor_expression_preview.R`. The controlled simulation’s preparation and probe steps precede the strong extension because the latter reuses the frozen base inputs and probe record. These commands describe the source workflow; they were not all executed as part of this package preparation.

The current main Figures 1–3 can be rendered from the frozen input tables in the same writable copy:

```sh
Rscript scripts/render_main_figures_1_2.R
Rscript scripts/render_main_figure_3.R
```

These scripts write `Figure_1`, `Figure_2`, and `Figure_3` PDF, SVG, and 600-dpi PNG files under `figures/`; refreshed panel tables under `outputs/figure_plotted_data/`; and QA/source maps under `outputs/figure_qa/` and `outputs/figure_source_manifests/`. They check the frozen numerical values and use `pdfinfo`, `pdftotext`, and `pdftoppm` from Poppler for PDF checks and previews.

The final figure renderer can be smoke-tested from the already bundled frozen tables without refitting models or rerunning simulations:

```sh
package_root="$(pwd)"
VS_OUTPUT_ROOT="$package_root" Rscript scripts/render_publication_figures.R
```

It writes `Figure_S15` through `Figure_S18` in PDF, PNG, and SVG formats under `figures/`.

## Software requirements

Use R 4.5.2 as the recorded execution version where possible. The scripts use the packages listed in `DEPENDENCIES.md`; the most demanding stages require Bioconductor packages MuSiC, SingleCellExperiment, S4Vectors, limma, and DESeq2. Figure generation requires ggplot2, patchwork, svglite, ragg, systemfonts, and an installed Arial font. Main-figure PDF checks require Poppler utilities on `PATH`; hash manifests call `shasum -a 256`.

`RUNTIME_VERSIONS.csv` records observed versions from the prior execution and a later host inspection; it is evidence, not an installation lock. No `renv.lock` is included, so a clean environment and full end-to-end portability have not been established.

## What this preparation did and did not verify

All 22 frozen input tables were copied byte-for-byte and hashed. Available generated R2/R3/R4/R5 tables and simulation prerequisites were added, with historical manifest paths rewritten to package-relative paths and checksums recomputed over the exact packaged bytes. Numerical tables were copied without value changes. The prior scenario-centers table has four column labels normalized from `r2_*` to `bulk_*`; its row order and values are unchanged and this is recorded in `metadata/SOURCE_MAP.csv`.

All 14 included R scripts passed syntax parsing. The supplementary renderer ran in a separate temporary work copy and produced all 12 S15–S18 PDF/PNG/SVG derivatives; PNG metadata recorded 600 dpi and the PDF font scan found embedded Arial-family fonts.

The two adapted main-figure renderers ran successfully in an isolated writable package copy. They produced all three one-page PDFs, SVGs, and 600-dpi PNGs with Arial-family PDF fonts. All ten generated panel plotting-data CSVs matched their included frozen snapshots byte-for-byte, and the built-in Fig. 1–3 numerical/value checks passed. The final dimensions were 173.9 × 175.0 mm (Figure 1), 173.9 × 159.8 mm (Figure 2), and 173.9 × 197.9 mm (Figure 3). This smoke test used Poppler's `pdfinfo`, `pdftotext`, and `pdftoppm` utilities; the R package `pdftools` is not required by these adapted scripts.

On 2026-10-19, the packaged bulk-model/effect-agreement script completed in an isolated writable copy, including its fixed 500-draw bootstrap. All 33 bulk, agreement, and S15/S16 result tables matched the frozen package tables exactly; path/hash-only metadata records were excluded from that comparison. This included 45,500 bootstrap metric rows and 35,500 resample-index rows.

The donor-reference script also completed in the isolated copy, fitting all six compartments. Its 118,716-row DESeq2 effects table and donor-count tables matched the frozen tables exactly. Common values in the other donor outputs also matched exactly. Differences were metadata-only: six tables contain an additional `cell_type` column, nine Schwann `lineage_id` values use `tumor` instead of the frozen `tumour` spelling, and one methods-description sentence was expanded. The donor run emitted DESeq2 line-search warnings but completed with zero raw-beta or apeglm convergence failures in its convergence output.

The controlled and strong simulation estimators and their validators were not rerun. A clean-environment, end-to-end run could further strengthen portability evidence; it is not presented here as a journal submission requirement. The limited model reruns and renderer smoke test do not establish that every stage is portable across environments.

The current package omits transient execution logs and superseded/debug attempts. It also does not include every intermediate record generated by the excluded R0/R1 accounting and document-build workflows. `INPUT_REQUIREMENTS.csv` lists remaining generated records by stage. The released package retains the tissue-only VIF gate failure, the 12/13 strong-validation warning, and the documented scope limits; these should be considered when interpreting the results.

## Citation and release metadata

This package is distributed as GitHub release v1.3.0 and archived through the repository's GitHub–Zenodo integration under concept DOI [10.5281/zenodo.21361686](https://doi.org/10.5281/zenodo.21361686). Cite the version-specific DOI shown on the Zenodo record for v1.3.0 when referring to this snapshot. The existing MIT license applies to code; it does not replace upstream data terms.
