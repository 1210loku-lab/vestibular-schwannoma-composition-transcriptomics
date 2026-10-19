#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
ROOT <- Sys.getenv("VS_INPUT_ROOT", unset = ".")
EXPORT_ROOT <- Sys.getenv("VS_OUTPUT_ROOT", unset = ".")
A <- file.path(EXPORT_ROOT, "outputs")
cts <- c("Schwann/肿瘤细胞", "巨噬/髓系", "T细胞", "成纤维", "内皮", "周细胞")
read_csv <- function(p) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
read_gz <- function(p) {
  con <- gzfile(p, open = "rt")
  on.exit(close(con), add = TRUE)
  utils::read.csv(con, stringsAsFactors = FALSE, check.names = FALSE)
}
write_csv <- function(x, p) utils::write.csv(x, p, row.names = FALSE, na = "", quote = TRUE)
sha256 <- function(p) {
  z <- system2("shasum", c("-a", "256", p), stdout = TRUE, stderr = TRUE)
  if (length(z) != 1L || !grepl("^[0-9a-fA-F]{64}\\s", z)) stop("SHA-256 failed: ", p)
  strsplit(z, "\\s+")[[1]][1]
}
checks <- list()
add <- function(name, pass, detail) checks[[length(checks) + 1L]] <<- data.frame(
  check = name, pass = isTRUE(pass), detail = as.character(detail), stringsAsFactors = FALSE
)

input_manifest <- read_csv(file.path(A, "simulation_FROZEN_INPUT_MANIFEST.csv"))
ih <- vapply(input_manifest$path, sha256, character(1))
add("frozen_input_sha256", all(ih == input_manifest$sha256), sprintf("%d/%d match", sum(ih == input_manifest$sha256), nrow(input_manifest)))
output_manifest <- read_csv(file.path(A, "simulation_OUTPUT_MANIFEST.csv"))
oh <- vapply(output_manifest$path, sha256, character(1))
add("runner_output_sha256", all(oh == output_manifest$sha256), sprintf("%d/%d match", sum(oh == output_manifest$sha256), nrow(output_manifest)))

truth <- read_csv(file.path(A, "simulation_gene_truth_and_marker_allocation.csv"))
profiles_tab <- read_gz(file.path(A, "simulation_baseline_linear_profiles.csv.gz"))
profiles <- as.matrix(profiles_tab[, -1, drop = FALSE])
storage.mode(profiles) <- "double"
rownames(profiles) <- profiles_tab[[1]]
classes <- table(truth$gene_class)
add("fixed_1500_gene_allocation", nrow(truth) == 1500L && !anyDuplicated(truth$gene) &&
      all(classes[c("null", "composition_only", "expression_only", "both")] == 100L) &&
      classes["reserved_marker"] == 300L && classes["non_target_background"] == 800L,
    paste(names(classes), as.integer(classes), sep = "=", collapse = "; "))
add("markers_disjoint_from_targets", !any(truth$reserved_marker &
      truth$gene_class %in% c("null", "composition_only", "expression_only", "both")),
    "reserved marker IDs do not overlap target evaluation classes")
pos <- vapply(c("expression_only", "both"), function(cl) sum(truth$gene_class == cl & truth$beta == 1.5), numeric(1))
neg <- vapply(c("expression_only", "both"), function(cl) sum(truth$gene_class == cl & truth$beta == -1.5), numeric(1))
add("fixed_signed_effects", all(pos == 50L) && all(neg == 50L),
    sprintf("positive counts=%s; negative counts=%s", paste(pos, collapse = "/"), paste(neg, collapse = "/")))
algebra <- read_csv(file.path(A, "simulation_noiseless_truth_identity_check.csv"))
max_algebra <- max(algebra$absolute_error)
add("noiseless_fixed_composition_identity", is.finite(max_algebra) && max_algebra < 1e-12,
    sprintf("max absolute error over %d genes=%.3g", nrow(algebra), max_algebra))
homogeneous <- profiles[truth$gene[truth$gene_class %in% c("null", "expression_only")], , drop = FALSE]
add("homogeneous_target_profiles", max(apply(homogeneous, 1, stats::sd)) < 1e-12,
    "null and expression-only profiles are constant across compartments")
add("equal_baseline_profile_totals", max(abs(colSums(profiles) - 10000)) < 1e-7,
    sprintf("maximum total error=%.3g", max(abs(colSums(profiles) - 10000))))

resource <- read_csv(file.path(A, "simulation_repeat_resource_metrics.csv"))
design <- read_csv(file.path(A, "simulation_design_and_gate_diagnostics.csv"))
fits <- read_gz(file.path(A, "simulation_gene_fit_outputs.csv.gz"))
warnings <- read_csv(file.path(A, "simulation_all_warnings_and_failures.csv"))
rep_metrics <- read_csv(file.path(A, "simulation_repeat_gene_class_metrics.csv"))
summary_metrics <- read_csv(file.path(A, "simulation_pilot_metric_summary.csv"))
fractions <- read_gz(file.path(A, "simulation_MuSiC_fraction_estimates_and_truth.csv.gz"))
frac_metrics <- read_csv(file.path(A, "simulation_fraction_recovery_metrics_by_repeat.csv"))
changes <- read_gz(file.path(A, "simulation_gene_effect_se_changes_vs_naive.csv.gz"))
scenario_counts <- table(resource$scenario)
add("20_repeats_per_scenario", nrow(resource) == 40L && length(scenario_counts) == 2L && all(scenario_counts == 20L),
    paste(names(scenario_counts), as.integer(scenario_counts), sep = "=", collapse = "; "))
add("three_routes_for_each_repeat", nrow(design) == 120L && all(table(design$route) == 40L),
    sprintf("design rows=%d; routes=%s", nrow(design), paste(table(design$route), collapse = "/")))
music_ok <- design$route == "estimated_music_composition" & design$music_status == "ok"
add("actual_music_success", sum(music_ok) == 40L && all(resource$music_status == "ok"), sprintf("MuSiC success=%d/40", sum(music_ok)))
add("all_limma_fits_complete", nrow(fits) == 180000L && all(table(fits$route) == 60000L) && all(fits$fit_status == "fit_ok"),
    sprintf("gene rows=%d; statuses=%s", nrow(fits), paste(names(table(fits$fit_status)), table(fits$fit_status), sep = "=", collapse = "; ")))
add("warnings_and_failures_retained", nrow(warnings) == 0L || all(nzchar(warnings$message)), sprintf("records=%d", nrow(warnings)))
add("runtime_and_heap_observed", all(is.finite(resource$runtime_seconds) & resource$runtime_seconds > 0) &&
      all(is.finite(resource$gc_peak_heap_mb_observed) & resource$gc_peak_heap_mb_observed > 0),
    sprintf("seconds %.3f..%.3f; R heap MB %.2f..%.2f", min(resource$runtime_seconds), max(resource$runtime_seconds),
            min(resource$gc_peak_heap_mb_observed), max(resource$gc_peak_heap_mb_observed)))
add("original_design_gate_columns", all(c("full_rank", "condition_number_unscaled", "group_vif", "all_original_gates_pass") %in% names(design)),
    sprintf("gate failures=%d/120", sum(!design$all_original_gates_pass, na.rm = TRUE)))
add("unique_gene_fit_keys", nrow(fits) == nrow(unique(fits[c("scenario", "replicate", "route", "gene")])),
    "one row per scenario × replicate × route × gene")
recalc_call <- !is.na(fits$fdr) & fits$fdr < 0.05 & abs(fits$logFC) > 1
recalc_direction <- fits$beta != 0 & sign(fits$logFC) != sign(fits$beta)
add("calls_and_directions_recomputed", all(recalc_call == fits$call_fdr05_abs_gt1) && all(recalc_direction == fits$direction_error),
    "saved calls and direction errors match independent recalculation")
add("metric_output_rows", nrow(rep_metrics) == 720L && nrow(summary_metrics) == 36L && nrow(changes) == 120000L,
    sprintf("repeat metrics=%d; summary=%d; gene effect/SE changes=%d", nrow(rep_metrics), nrow(summary_metrics), nrow(changes)))

# Recompute binary-rate MCSE from the repeat-level rates.
mcse_ok <- TRUE
for (i in seq_len(nrow(summary_metrics))) {
  s <- summary_metrics[i, ]
  x <- rep_metrics$threshold_call_rate_all[rep_metrics$scenario == s$scenario & rep_metrics$route == s$route & rep_metrics$gene_class == s$gene_class]
  x <- x[is.finite(x)]
  expected <- if (length(x) > 1L) stats::sd(x) / sqrt(length(x)) else NA_real_
  if (is.finite(expected) && (!is.finite(s$call_rate_mcse) || abs(expected - s$call_rate_mcse) > 1e-10)) mcse_ok <- FALSE
}
add("binary_rate_MCSE_recomputed", mcse_ok, "SD of repeat-level rates divided by sqrt(finite repeats)")

truth_fraction <- read_gz(file.path(A, "simulation_mixture_fraction_truth.csv.gz"))
centers <- read_csv(file.path(A, "simulation_scenario_centers_and_bulk_contrast.csv"))
sep_centers <- centers[centers$scenario == "separation", , drop = FALSE]
bulk_comp <- read_csv(file.path(A, "bulk_composition_distribution.csv"))
bulk_music <- bulk_comp[bulk_comp$reference == "MuSiC" & bulk_comp$feature %in% c("Schwann", "Myeloid", "T_cell", "Fibroblast", "Endothelial", "Pericyte") & bulk_comp$group %in% c("Control", "Tumor"), , drop = FALSE]
bulk_control <- bulk_music[bulk_music$group == "Control", , drop = FALSE]
bulk_tumor <- bulk_music[bulk_music$group == "Tumor", , drop = FALSE]
bulk_control <- bulk_control$mean[match(c("Schwann", "Myeloid", "T_cell", "Fibroblast", "Endothelial", "Pericyte"), bulk_control$feature)]
bulk_tumor <- bulk_tumor$mean[match(c("Schwann", "Myeloid", "T_cell", "Fibroblast", "Endothelial", "Pericyte"), bulk_tumor$feature)]
bulk_control <- bulk_control / sum(bulk_control)
bulk_tumor <- bulk_tumor / sum(bulk_tumor)
center_match <- nrow(sep_centers) == length(cts) &&
  max(abs(sep_centers$control_center - bulk_control)) < 1e-12 &&
  max(abs(sep_centers$tumor_center - bulk_tumor)) < 1e-12
center_difference <- if (nrow(sep_centers) == length(cts)) {
  max(abs(sep_centers$center_tumor_minus_control))
} else NA_real_
truth_difference <- if (nrow(sep_centers) == length(cts)) {
  max(abs(sep_centers$realized_mean_cell_fraction_tumor_minus_control))
} else NA_real_
add("corrected_separation_centers_match_bulk", center_match && is.finite(center_difference) && center_difference > 0.10,
    sprintf("bulk-model normalized centers exact=%s; maximum component contrast=%.6f", center_match, center_difference))
add("separation_retained_in_realized_cell_truth", is.finite(truth_difference) && truth_difference > 0.10,
    sprintf("maximum absolute realized cell-fraction group contrast=%.6f", truth_difference))
sum_cell <- aggregate(cell_fraction ~ scenario + replicate + sample, truth_fraction, sum)
sum_rna <- aggregate(rna_contribution_fraction ~ scenario + replicate + sample, truth_fraction, sum)
sum_music <- aggregate(music_estimated_fraction ~ scenario + replicate + sample, fractions, sum)
add("separate_cell_and_RNA_truth", max(abs(sum_cell$cell_fraction - 1)) < 1e-10 && max(abs(sum_rna$rna_contribution_fraction - 1)) < 1e-10,
    sprintf("cell sum error=%.3g; RNA sum error=%.3g", max(abs(sum_cell$cell_fraction - 1)), max(abs(sum_rna$rna_contribution_fraction - 1))))
add("fraction_estimates_complete", nrow(fractions) == 9600L && all(is.finite(fractions$music_estimated_fraction)) &&
      max(abs(sum_music$music_estimated_fraction - 1)) < 1e-8 && nrow(frac_metrics) == 40L,
    sprintf("estimate rows=%d; max sum error=%.3g", nrow(fractions), max(abs(sum_music$music_estimated_fraction - 1))))

# Independent direct MuSiC call on one fixed replicate verifies saved fractions.
if (!requireNamespace("MuSiC", quietly = TRUE) || !requireNamespace("SingleCellExperiment", quietly = TRUE) ||
    !requireNamespace("S4Vectors", quietly = TRUE) || !requireNamespace("limma", quietly = TRUE)) stop("Validation package missing")
suppressPackageStartupMessages({
  library(MuSiC)
  library(SingleCellExperiment)
  library(S4Vectors)
  library(limma)
})
ref_tab <- read_gz(file.path(A, "simulation_synthetic_reference_counts.csv.gz"))
ref_counts <- as.matrix(ref_tab[, -1, drop = FALSE]); storage.mode(ref_counts) <- "double"
rownames(ref_counts) <- ref_tab[[1]]
ref_meta <- read_csv(file.path(A, "simulation_synthetic_reference_cell_metadata.csv"))
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = ref_counts),
  colData = S4Vectors::DataFrame(ref_meta[, c("sample", "cell_type"), drop = FALSE])
)
mix_tab <- read_gz(file.path(A, "simulation_mixture_log2_expression.csv.gz"))
mix <- as.matrix(mix_tab[, -1, drop = FALSE]); storage.mode(mix) <- "double"
rownames(mix) <- mix_tab[[1]]
sample_meta <- read_csv(file.path(A, "simulation_mixture_sample_metadata.csv"))
verify_meta <- sample_meta[sample_meta$scenario == "overlap" & sample_meta$replicate == 1L, , drop = FALSE]
ids <- verify_meta$sample
markers <- truth$gene[truth$reserved_marker]
music_warnings <- character()
verify <- withCallingHandlers(MuSiC::music_prop(
  bulk.mtx = 2^mix[, ids, drop = FALSE], sc.sce = sce, markers = markers,
  clusters = "cell_type", samples = "sample", verbose = FALSE,
  iter.max = 1000, nu = 1e-4, eps = 0.01, centered = FALSE, normalize = TRUE
), warning = function(w) {
  music_warnings <<- c(music_warnings, conditionMessage(w))
  invokeRestart("muffleWarning")
})
music_direct <- as.matrix(verify$Est.prop.weighted)[ids, cts, drop = FALSE]
saved_fraction <- fractions[fractions$scenario == "overlap" & fractions$replicate == 1L, , drop = FALSE]
music_saved <- matrix(NA_real_, length(ids), length(cts), dimnames = list(ids, cts))
for (sid in ids) {
  z <- saved_fraction[saved_fraction$sample == sid, , drop = FALSE]
  music_saved[sid, ] <- z$music_estimated_fraction[match(cts, z$cell_type)]
}
music_diff <- max(abs(music_direct - music_saved))
add("independent_actual_MuSiC_fraction_check", is.finite(music_diff) && music_diff < 1e-8 && !length(music_warnings),
    sprintf("overlap replicate 1 max difference=%.3g; warnings=%d", music_diff, length(music_warnings)))

# Independently recompute one naive coefficient and unmoderated SE using stats::lm.
effect_gene <- truth$gene[truth$gene_class == "expression_only" & truth$beta == 1.5][1]
y <- mix[effect_gene, ids]
group <- factor(verify_meta$group, levels = c("Control", "Tumor"))
design_one <- stats::model.matrix(~ group)
limma_one <- limma::lmFit(matrix(y, nrow = 1L, dimnames = list(effect_gene, ids)), design_one)
limma_beta <- unname(limma_one$coefficients[1, "groupTumor"])
limma_se <- unname(limma_one$stdev.unscaled[1, "groupTumor"] * limma_one$sigma[1])
ols <- stats::lm(y ~ group)
ols_beta <- unname(stats::coef(ols)["groupTumor"])
ols_se <- unname(summary(ols)$coefficients["groupTumor", "Std. Error"])
saved_gene <- fits[fits$scenario == "overlap" & fits$replicate == 1L & fits$route == "naive" & fits$gene == effect_gene, , drop = FALSE]
fit_diff <- max(abs(c(limma_beta - ols_beta, limma_se - ols_se,
                      limma_beta - saved_gene$logFC, limma_se - saved_gene$raw_se)))
add("independent_limma_coefficient_and_SE_check", nrow(saved_gene) == 1L && fit_diff < 1e-9,
    sprintf("gene=%s; max lmFit/OLS/saved discrepancy=%.3g", effect_gene, fit_diff))

validation <- do.call(rbind, checks)
write_csv(validation, file.path(A, "simulation_VALIDATION.csv"))
cat("simulation independent output validation\n")
cat(paste(sprintf("%s: %s — %s", validation$check, ifelse(validation$pass, "PASS", "FAIL"), validation$detail), collapse = "\n"), "\n")
cat("checks=", nrow(validation), "; passed=", sum(validation$pass), "; failed=", sum(!validation$pass), "\n", sep = "")
if (any(!validation$pass)) quit(save = "no", status = 1L)
