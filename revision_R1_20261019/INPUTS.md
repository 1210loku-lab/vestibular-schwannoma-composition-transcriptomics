# Input and frozen-output map

`VS_INPUT_ROOT` must point to this package's `inputs/` directory. The 13 processed model/simulation inputs and nine additional main-figure source tables preserve their source-relative `data/processed/` and `results/` paths below `inputs/`. Together, 22 frozen input tables are listed with byte counts and SHA-256 values in `metadata/INPUT_SHA256.csv`; Figure 1 reuses the already-bundled `deg_composition_adjusted.csv` table.

The included model and simulation scripts write into `VS_OUTPUT_ROOT/outputs/`. The main-figure renderers write vector and raster figures to `VS_OUTPUT_ROOT/figures/`, refreshed panel plotting tables to `VS_OUTPUT_ROOT/outputs/figure_plotted_data/`, QA outputs to `VS_OUTPUT_ROOT/outputs/figure_qa/`, and package-relative source maps to `VS_OUTPUT_ROOT/outputs/figure_source_manifests/`. Frozen plotted-data snapshots are included for comparison with the current Fig. 1–3 renderers.

This package also carries the archived frozen prerequisites and fit/metric tables needed to inspect prior results. The simulation manifests in `outputs/` use package-relative paths so their frozen hashes can be checked from the package root. These packaged manifests are integrity records for bundled bytes; they are not new pre-fit approvals or a new statistical validation. No Fig. 4 or provisional UMAP source is included.

See `INPUT_REQUIREMENTS.csv` for each script input, generated prerequisite, bundled result, or intentionally omitted log. See `metadata/SOURCE_MAP.csv` for source provenance and transformations.
