#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
ROOT <- Sys.getenv("VS_INPUT_ROOT", unset = ".")
EXPORT_ROOT <- Sys.getenv("VS_OUTPUT_ROOT", unset = ".")
A <- file.path(EXPORT_ROOT, "outputs")
SCRIPT <- file.path(EXPORT_ROOT, "scripts", "probe_all_gene_strong_simulation.R")
cell_types <- c("Schwann/肿瘤细胞", "巨噬/髓系", "T细胞", "成纤维", "内皮", "周细胞")
features <- c("Schwann", "Myeloid", "T_cell", "Fibroblast", "Endothelial", "Pericyte")
read_csv <- function(p) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
read_gz <- function(p) {
  con <- gzfile(p, open = "rt")
  on.exit(close(con), add = TRUE)
  utils::read.csv(con, stringsAsFactors = FALSE, check.names = FALSE)
}
write_csv <- function(x, p) utils::write.csv(x, p, row.names = FALSE, na = "", quote = TRUE)
sha256 <- function(p) {
  z <- system2("shasum", c("-a", "256", p), stdout = TRUE, stderr = TRUE)
  if (length(z) != 1L || !grepl("^[0-9a-fA-F]{64}\\s", z)) stop("SHA-256 failed for ", p)
  strsplit(z, "\\s+")[[1]][1]
}

plan <- file.path(A, "simulation_strong_ALLGENE_PROBE_PLAN.md")
protocol <- file.path(A, "simulation_strong_PROTOCOL.json")
input_manifest_path <- file.path(A, "simulation_strong_FROZEN_INPUT_MANIFEST.csv")
input_manifest <- read_csv(input_manifest_path)
if (!all(file.exists(input_manifest$path))) stop("Frozen primary-pilot inputs are missing")
input_hashes <- vapply(input_manifest$path, sha256, character(1))
if (any(input_hashes != input_manifest$sha256)) stop("Primary-pilot frozen input hash mismatch")
probe_input_paths <- c(plan, protocol, input_manifest_path,
  file.path(A, "simulation_synthetic_reference_counts.csv.gz"),
  file.path(A, "simulation_synthetic_reference_cell_metadata.csv"),
  file.path(A, "simulation_strong_mixture_log2_expression.csv.gz"),
  file.path(A, "simulation_strong_sample_metadata.csv"),
  file.path(A, "simulation_strong_fraction_truth.csv.gz"), SCRIPT)
probe_input_manifest <- data.frame(
  path = normalizePath(probe_input_paths, winslash = "/", mustWork = TRUE),
  bytes = file.info(probe_input_paths)$size,
  sha256 = vapply(probe_input_paths, sha256, character(1)),
  stringsAsFactors = FALSE
)
write_csv(probe_input_manifest, file.path(A, "simulation_strong_ALLGENE_PROBE_INPUT_MANIFEST.csv"))

if (!requireNamespace("MuSiC", quietly = TRUE) || !requireNamespace("SingleCellExperiment", quietly = TRUE) ||
    !requireNamespace("S4Vectors", quietly = TRUE)) stop("All-gene probe packages are missing")
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(S4Vectors)
})
ref_tab <- read_gz(file.path(A, "simulation_synthetic_reference_counts.csv.gz"))
ref_counts <- as.matrix(ref_tab[, -1, drop = FALSE])
storage.mode(ref_counts) <- "double"
rownames(ref_counts) <- ref_tab[[1]]
ref_meta <- read_csv(file.path(A, "simulation_synthetic_reference_cell_metadata.csv"))
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = ref_counts),
  colData = S4Vectors::DataFrame(ref_meta[, c("sample", "cell_type"), drop = FALSE])
)
mix_tab <- read_gz(file.path(A, "simulation_strong_mixture_log2_expression.csv.gz"))
bulk <- as.matrix(mix_tab[, -1, drop = FALSE])
storage.mode(bulk) <- "double"
rownames(bulk) <- mix_tab[[1]]
sample_meta <- read_csv(file.path(A, "simulation_strong_sample_metadata.csv"))
truth <- read_gz(file.path(A, "simulation_strong_fraction_truth.csv.gz"))

probe_rows <- list()
estimate_rows <- list()
log_lines <- c(paste("simulation strong all-common-gene MuSiC feasibility probe; package", as.character(utils::packageVersion("MuSiC"))),
  "Fixed selection: replicate 1, one call per frozen scenario, markers=NULL, normalize=TRUE, iter.max=1000, nu=1e-4, eps=0.01, centered=FALSE.")
for (scenario_name in c("strong")) {
  meta <- sample_meta[sample_meta$scenario == scenario_name & sample_meta$replicate == 1L, , drop = FALSE]
  ids <- meta$sample
  warning_text <- character()
  started <- proc.time()[[3]]
  res <- tryCatch(withCallingHandlers(
    MuSiC::music_prop(
      bulk.mtx = 2^bulk[, ids, drop = FALSE], sc.sce = sce,
      markers = NULL, clusters = "cell_type", samples = "sample", verbose = FALSE,
      iter.max = 1000, nu = 1e-4, eps = 0.01, centered = FALSE, normalize = TRUE
    ),
    warning = function(w) {
      warning_text <<- c(warning_text, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  ), error = function(e) e)
  elapsed <- proc.time()[[3]] - started
  gc_observed <- gc()
  heap_mb <- sum(gc_observed[, 7])
  status <- if (inherits(res, "error")) "error" else "ok"
  error_text <- if (inherits(res, "error")) conditionMessage(res) else ""
  log_lines <- c(log_lines, sprintf("%s replicate 1: status=%s; runtime=%.3f sec; warnings=%d; error=%s",
                                    scenario_name, status, elapsed, length(warning_text), error_text))
  if (length(warning_text)) log_lines <- c(log_lines, paste0("warning: ", warning_text))
  if (status == "ok") {
    est <- as.matrix(res$Est.prop.weighted)[ids, cell_types, drop = FALSE]
    colnames(est) <- features
    if (any(!is.finite(est))) {
      status <- "invalid_fraction_values"
      error_text <- "MuSiC returned non-finite all-gene fractions"
    } else {
      truth_sub <- truth[truth$scenario == scenario_name & truth$replicate == 1L, , drop = FALSE]
      truth_matrix <- matrix(NA_real_, nrow = length(ids), ncol = length(features), dimnames = list(ids, features))
      rna_matrix <- truth_matrix
      for (sid in ids) {
        z <- truth_sub[truth_sub$sample == sid, , drop = FALSE]
        truth_matrix[sid, ] <- z$cell_fraction[match(cell_types, z$cell_type)]
        rna_matrix[sid, ] <- z$rna_contribution_fraction[match(cell_types, z$cell_type)]
        for (j in seq_along(features)) {
          estimate_rows[[length(estimate_rows) + 1L]] <- data.frame(
            scenario = scenario_name, replicate = 1L, sample = sid,
            group = meta$group[match(sid, meta$sample)], cell_type = features[j],
            allgene_estimate = est[sid, j], cell_fraction_truth = truth_matrix[sid, j],
            rna_contribution_truth = rna_matrix[sid, j],
            abs_error_vs_cell_fraction = abs(est[sid, j] - truth_matrix[sid, j]),
            abs_error_vs_rna_contribution = abs(est[sid, j] - rna_matrix[sid, j]),
            stringsAsFactors = FALSE
          )
        }
      }
      probe_rows[[length(probe_rows) + 1L]] <- data.frame(
        scenario = scenario_name, replicate = 1L, n_samples = length(ids),
        marker_argument = "NULL_all_nonzero_common_genes", status = status,
        runtime_seconds = elapsed, gc_peak_heap_mb_observed = heap_mb,
        warning_count = length(warning_text), error = error_text,
        max_fraction_sum_error = max(abs(rowSums(est) - 1)),
        mean_absolute_error_vs_cell_fraction = mean(abs(est - truth_matrix)),
        rmse_vs_cell_fraction = sqrt(mean((est - truth_matrix)^2)),
        mean_absolute_error_vs_rna_contribution = mean(abs(est - rna_matrix)),
        stringsAsFactors = FALSE
      )
    }
  }
  if (status != "ok") {
    probe_rows[[length(probe_rows) + 1L]] <- data.frame(
      scenario = scenario_name, replicate = 1L, n_samples = length(ids),
      marker_argument = "NULL_all_nonzero_common_genes", status = status,
      runtime_seconds = elapsed, gc_peak_heap_mb_observed = heap_mb,
      warning_count = length(warning_text), error = error_text,
      max_fraction_sum_error = NA_real_, mean_absolute_error_vs_cell_fraction = NA_real_,
      rmse_vs_cell_fraction = NA_real_, mean_absolute_error_vs_rna_contribution = NA_real_,
      stringsAsFactors = FALSE
    )
  }
}
probe_table <- do.call(rbind, probe_rows)
estimate_table <- if (length(estimate_rows)) do.call(rbind, estimate_rows) else data.frame()
write_csv(probe_table, file.path(A, "simulation_strong_ALLGENE_PROBE_SUMMARY.csv"))
write_csv(estimate_table, file.path(A, "simulation_strong_ALLGENE_PROBE_FRACTIONS.csv"))
writeLines(log_lines, file.path(A, "simulation_strong_ALLGENE_PROBE.log"))
output_paths <- c(SCRIPT, file.path(A, "simulation_strong_ALLGENE_PROBE_PLAN.md"),
  file.path(A, "simulation_strong_ALLGENE_PROBE_INPUT_MANIFEST.csv"), file.path(A, "simulation_strong_ALLGENE_PROBE_SUMMARY.csv"),
  file.path(A, "simulation_strong_ALLGENE_PROBE_FRACTIONS.csv"), file.path(A, "simulation_strong_ALLGENE_PROBE.log"))
out_manifest <- data.frame(
  path = normalizePath(output_paths, winslash = "/", mustWork = TRUE),
  bytes = file.info(output_paths)$size,
  sha256 = vapply(output_paths, sha256, character(1)),
  stringsAsFactors = FALSE
)
write_csv(out_manifest, file.path(A, "simulation_strong_ALLGENE_PROBE_OUTPUT_MANIFEST.csv"))
cat("simulation strong all-common-gene probe complete; calls=", nrow(probe_table),
    "; successes=", sum(probe_table$status == "ok"), "; runtime seconds=",
    paste(round(probe_table$runtime_seconds, 3), collapse = "/"), "\n", sep = "")
if (any(probe_table$status != "ok")) quit(save = "no", status = 1L)
