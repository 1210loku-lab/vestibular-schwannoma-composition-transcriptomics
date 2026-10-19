# R dependencies

No dependency lockfile is included. The following is a static summary of packages loaded or required by the bundled scripts; versions recorded during prior execution are in `RUNTIME_VERSIONS.csv`.

| Script group | Packages |
|---|---|
| Bulk model and effect agreement | limma, data.table, ggplot2, patchwork, openssl, svglite, ragg |
| Donor-reference analysis and preview | DESeq2, ggplot2 |
| Controlled/strong simulation, probes, validation | MuSiC, SingleCellExperiment, S4Vectors, limma, data.table, jsonlite |
| Main Figures 1–3 | ggplot2, patchwork, svglite, ragg, systemfonts, scales |
| Final supplementary figures (S15–S18) | ggplot2, patchwork, svglite, ragg |

Base R namespaces used include `stats`, `utils`, `grDevices`, and `grid`. The figure scripts request Arial; no font file is bundled. Several manifest writers call `shasum -a 256`, so a compatible `shasum` command must be on `PATH`.

The main-figure renderers use Poppler command-line utilities `pdfinfo`, `pdftotext`, and `pdftoppm` for PDF page checks, text checks, and QA previews. These executables must be on `PATH`; their version is recorded in `RUNTIME_VERSIONS.csv` as host evidence.
