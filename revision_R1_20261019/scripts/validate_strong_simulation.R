#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE, warn = 1)
ROOT <- Sys.getenv("VS_INPUT_ROOT", unset = ".")
EXPORT_ROOT <- Sys.getenv("VS_OUTPUT_ROOT", unset = ".")
A <- file.path(EXPORT_ROOT, "outputs")
cell_types <- c("Schwann/肿瘤细胞", "巨噬/髓系", "T细胞", "成纤维", "内皮", "周细胞")
features <- c("Schwann", "Myeloid", "T_cell", "Fibroblast", "Endothelial", "Pericyte")
read_csv <- function(p) utils::read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
read_gz <- function(p) {
  con <- gzfile(p, open = "rt"); on.exit(close(con), add = TRUE)
  utils::read.csv(con, check.names = FALSE, stringsAsFactors = FALSE)
}
sha256 <- function(p) {
  z <- system2("shasum", c("-a", "256", p), stdout = TRUE, stderr = TRUE)
  if (length(z) != 1L || !grepl("^[0-9a-fA-F]{64}\\s", z)) stop("SHA-256 failed for ", p)
  strsplit(z, "\\s+")[[1]][1]
}
rec <- function(name, pass, detail) data.frame(check = name, pass = isTRUE(pass), detail = detail, stringsAsFactors = FALSE)
checks <- list()

review <- jsonlite::fromJSON(file.path(EXPORT_ROOT, "outputs", "protocol_integrity_lock.json"), simplifyVector = FALSE)
protocol_hash <- sha256(file.path(A, "simulation_strong_PROTOCOL.json"))
checks[[length(checks) + 1L]] <- rec("protocol_integrity_lock_matches", identical(protocol_hash, review$protocol_sha256), protocol_hash)
frozen <- read_csv(file.path(A, "simulation_strong_FROZEN_INPUT_MANIFEST.csv"))
frozen_hash <- vapply(frozen$path, sha256, character(1))
checks[[length(checks) + 1L]] <- rec("frozen_input_hashes_match", all(frozen_hash == frozen$sha256),
  sprintf("%d/%d frozen inputs match", sum(frozen_hash == frozen$sha256), nrow(frozen)))

gene_truth <- read_csv(file.path(A, "simulation_gene_truth_and_marker_allocation.csv"))
profile_tab <- read_gz(file.path(A, "simulation_baseline_linear_profiles.csv.gz"))
profile <- as.matrix(profile_tab[, -1, drop = FALSE]); storage.mode(profile) <- "double"; rownames(profile) <- profile_tab[[1]]
beta <- setNames(gene_truth$beta, gene_truth$gene)
mix_tab <- read_gz(file.path(A, "simulation_strong_mixture_log2_expression.csv.gz"))
mix <- as.matrix(mix_tab[, -1, drop = FALSE]); storage.mode(mix) <- "double"; rownames(mix) <- mix_tab[[1]]
meta <- read_csv(file.path(A, "simulation_strong_sample_metadata.csv"))
truth <- read_gz(file.path(A, "simulation_strong_fraction_truth.csv.gz"))
old_meta <- read_csv(file.path(A, "simulation_mixture_sample_metadata.csv"))
seeds <- sort(unique(as.integer(meta$seed)))
seed_ok <- length(seeds) == 20L && identical(seeds, 2026101801:2026101820) &&
  !length(intersect(seeds, as.integer(unique(old_meta$seed))))
group_counts <- table(meta$replicate, meta$group)
sample_ok <- nrow(meta) == 800L && ncol(mix) == 800L && nrow(truth) == 4800L &&
  identical(colnames(mix), meta$sample) && !anyDuplicated(meta$sample) &&
  all(group_counts[, "Control"] == 9L) && all(group_counts[, "Tumor"] == 31L)
checks[[length(checks) + 1L]] <- rec("20_fixed_new_seeds_and_800_sample_scope", seed_ok && sample_ok,
  sprintf("samples=%d; truth rows=%d; seeds=%d-%d; old-seed overlap=%d", nrow(meta), nrow(truth), min(seeds), max(seeds), length(intersect(seeds, as.integer(unique(old_meta$seed))))))
cell_sums <- aggregate(cell_fraction ~ sample, truth, sum)$cell_fraction
rna_sums <- aggregate(rna_contribution_fraction ~ sample, truth, sum)$rna_contribution_fraction
sum_error <- max(abs(cell_sums - 1), abs(rna_sums - 1))
checks[[length(checks) + 1L]] <- rec("cell_and_RNA_truth_separately_sum_to_one", sum_error < 1e-10, sprintf("max sum error=%.3g", sum_error))

centers <- read_csv(file.path(A, "simulation_strong_centers_and_truth_contrast.csv"))
base_centers <- read_csv(file.path(A, "simulation_scenario_centers_and_bulk_contrast.csv"))
moderate <- base_centers[base_centers$scenario == "separation", , drop = FALSE]
bulk_composition <- read_csv(file.path(A, "bulk_composition_distribution.csv"))
bulk_composition <- bulk_composition[bulk_composition$reference == "MuSiC" & bulk_composition$feature %in% features & bulk_composition$group %in% c("Control", "Tumor"), , drop = FALSE]
get_bulk <- function(group) {
  d <- bulk_composition[bulk_composition$group == group, , drop = FALSE]
  x <- d$mean[match(features, d$feature)]; x / sum(x)
}
center_err <- max(abs(centers$control_center - moderate$control_center), abs(centers$tumor_center - moderate$tumor_center),
                  abs(centers$control_center - get_bulk("Control")), abs(centers$tumor_center - get_bulk("Tumor")))
center_sep <- max(abs(centers$center_tumor_minus_control))
realized_sep <- max(abs(centers$realized_mean_cell_fraction_tumor_minus_control))
checks[[length(checks) + 1L]] <- rec("exact_bulk_centers_and_distinct_realized_truth", center_err < 1e-12 && center_sep > .1 && realized_sep > .1,
  sprintf("max center mismatch=%.3g; bulk-model center contrast=%.6f; realized mean cell-fraction contrast=%.6f", center_err, center_sep, realized_sep))
p0 <- truth[truth$sample == meta$sample[1], , drop = FALSE]
p0 <- p0$cell_fraction[match(features, p0$feature)]
linear0 <- drop(profile %*% p0)
identity_error <- max(abs(log2((linear0 * 2^unname(beta[rownames(profile)])) / linear0) - unname(beta[rownames(profile)])))
identity_saved <- read_csv(file.path(A, "simulation_strong_noiseless_truth_identity_check.csv"))
identity_saved_error <- max(abs(identity_saved$absolute_error))
checks[[length(checks) + 1L]] <- rec("noiseless_same_composition_log2_contrast_equals_beta",
  identity_error < 1e-12 && identity_saved_error < 1e-12,
  sprintf("independent max error=%.3g; saved max error=%.3g", identity_error, identity_saved_error))

out_manifest <- read_csv(file.path(A, "simulation_strong_OUTPUT_MANIFEST.csv"))
out_hash <- vapply(out_manifest$path, sha256, character(1))
checks[[length(checks) + 1L]] <- rec("output_manifest_hashes_match", all(out_hash == out_manifest$sha256),
  sprintf("%d/%d outputs match", sum(out_hash == out_manifest$sha256), nrow(out_manifest)))
fits <- read_gz(file.path(A, "simulation_strong_gene_fit_outputs.csv.gz"))
diag <- read_csv(file.path(A, "simulation_strong_design_and_gate_diagnostics.csv"))
warnings <- read_csv(file.path(A, "simulation_strong_all_warnings_and_failures.csv"))
routes <- c("naive", "oracle_true_composition", "estimated_music_composition")
fit_ok <- nrow(fits) == 90000L && nrow(diag) == 60L && all(table(fits$replicate, fits$route) == 1500L) &&
  setequal(unique(fits$replicate), 1:20) && setequal(unique(fits$route), routes) && all(diag$n_samples == 40L)
checks[[length(checks) + 1L]] <- rec("20x3_route_fit_rows_and_all_diagnostics_retained", fit_ok,
  sprintf("gene rows=%d; design rows=%d; warning/failure records=%d", nrow(fits), nrow(diag), nrow(warnings)))
vif <- diag$group_vif[diag$route == "estimated_music_composition"]
condition <- diag$condition_number_unscaled[diag$route == "estimated_music_composition"]
gate_fail <- sum(!diag$all_original_gates_pass, na.rm = TRUE)
checks[[length(checks) + 1L]] <- rec("original_design_gates_reported_without_filtering", length(vif) == 20L &&
  all(is.finite(vif)) && all(vif <= 10) && all(is.finite(condition)) && all(condition <= 1000),
  sprintf("all-route gate failures=%d; estimated VIF min/median/max=%.3f/%.3f/%.3f; condition max=%.3f", gate_fail, min(vif), median(vif), max(vif), max(condition)))

metrics <- read_csv(file.path(A, "simulation_strong_repeat_gene_class_metrics.csv"))
summary <- read_csv(file.path(A, "simulation_strong_pilot_metric_summary.csv"))
mcse_ok <- nrow(metrics) == 360L && nrow(summary) == 18L && all(summary$attempted_repeats == 20L) && all(summary$scenario == "strong")
mcse_diff <- 0
if (mcse_ok) {
  for (i in seq_len(nrow(summary))) {
    x <- metrics[metrics$route == summary$route[i] & metrics$gene_class == summary$gene_class[i], , drop = FALSE]
    for (field in c(false = "false_call_rate_beta_zero", sensitivity = "sensitivity_nonzero_beta")) {
      vals <- x[[field]]
      expected <- if (sum(is.finite(vals)) > 1L) stats::sd(vals, na.rm = TRUE) / sqrt(sum(is.finite(vals))) else NA_real_
      actual_name <- if (field == "false") "beta_zero_false_call_rate_mcse" else "nonzero_beta_sensitivity_mcse"
      actual <- summary[[actual_name]][i]
      d <- if (is.na(expected) && is.na(actual)) 0 else abs(expected - actual)
      mcse_diff <- max(mcse_diff, d, na.rm = TRUE)
    }
  }
}
checks[[length(checks) + 1L]] <- rec("target_metrics_include_independently_recomputed_20_repeat_MCSE", mcse_ok && mcse_diff < 1e-10,
  sprintf("repeat/class/route rows=%d; summaries=%d; max MCSE difference=%.3g", nrow(metrics), nrow(summary), mcse_diff))

# Exact-input MuSiC replay checks the saved estimated fractions for one frozen repeat.
suppressPackageStartupMessages({library(SingleCellExperiment); library(S4Vectors)})
ref_tab <- read_gz(file.path(A, "simulation_synthetic_reference_counts.csv.gz"))
ref_counts <- as.matrix(ref_tab[, -1, drop = FALSE]); storage.mode(ref_counts) <- "double"; rownames(ref_counts) <- ref_tab[[1]]
ref_meta <- read_csv(file.path(A, "simulation_synthetic_reference_cell_metadata.csv"))
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = ref_counts),
  colData = S4Vectors::DataFrame(ref_meta[, c("sample", "cell_type"), drop = FALSE]))
meta1 <- meta[meta$replicate == 1L, , drop = FALSE]; ids <- meta1$sample
markers <- gene_truth$gene[gene_truth$reserved_marker]
started <- proc.time()[[3]]
replay <- MuSiC::music_prop(bulk.mtx = 2^mix[, ids, drop = FALSE], sc.sce = sce, markers = markers,
  clusters = "cell_type", samples = "sample", verbose = FALSE, iter.max = 1000,
  nu = 1e-4, eps = 0.01, centered = FALSE, normalize = TRUE)
replay_seconds <- proc.time()[[3]] - started
gc_observed <- gc(); replay_heap <- sum(gc_observed[, 7])
est <- as.matrix(replay$Est.prop.weighted)[ids, cell_types, drop = FALSE]
est[est < 0] <- 0; est <- est / rowSums(est)
saved_fraction <- read_gz(file.path(A, "simulation_strong_MuSiC_fraction_estimates_and_truth.csv.gz"))
saved1 <- saved_fraction[saved_fraction$replicate == 1L, , drop = FALSE]
saved_est <- matrix(saved1$music_estimated_fraction[match(paste(rep(ids, each = length(cell_types)), rep(cell_types, times = length(ids))), paste(saved1$sample, saved1$cell_type))],
  nrow = length(ids), ncol = length(cell_types), byrow = TRUE, dimnames = list(ids, cell_types))
music_diff <- max(abs(est - saved_est))
checks[[length(checks) + 1L]] <- rec("independent_actual_MuSiC_replay_matches_repeat1", music_diff < 1e-8,
  sprintf("max fraction difference=%.3g; replay runtime=%.3f sec; observed R GC/vector heap=%.1f MB", music_diff, replay_seconds, replay_heap))

# Independent limma refit checks one coefficient and its unmoderated standard error.
z <- log(pmax(est, 1e-6)); clr <- z - rowMeans(z)
pc <- stats::prcomp(clr, center = TRUE, scale. = FALSE)$x[, 1:5, drop = FALSE]; colnames(pc) <- paste0("PC", 1:5)
dat <- data.frame(group = factor(meta1$group, levels = c("Control", "Tumor")), pc, check.names = FALSE)
design <- model.matrix(~ ., data = dat)
gene0 <- gene_truth$gene[which(gene_truth$gene_class == "both")[1]]
indfit <- limma::lmFit(mix[gene0, ids, drop = FALSE], design)
coef_i <- match("groupTumor", colnames(indfit$coefficients))
ind_beta <- as.numeric(indfit$coefficients[1, coef_i])
ind_se <- as.numeric(indfit$stdev.unscaled[1, coef_i] * indfit$sigma[1])
saved_gene <- fits[fits$replicate == 1L & fits$route == "estimated_music_composition" & fits$gene == gene0, , drop = FALSE]
beta_diff <- abs(ind_beta - saved_gene$logFC); se_diff <- abs(ind_se - saved_gene$raw_se)
checks[[length(checks) + 1L]] <- rec("independent_limma_beta_and_raw_SE_match", nrow(saved_gene) == 1L && beta_diff < 1e-8 && se_diff < 1e-8,
  sprintf("gene=%s; beta diff=%.3g; raw-SE diff=%.3g; beta=%.8f; known beta=%.3f", gene0, beta_diff, se_diff, ind_beta, gene_truth$beta[match(gene0, gene_truth$gene)]))

base_probe <- read_csv(file.path(A, "simulation_ALLGENE_PROBE_SUMMARY.csv"))
strong_probe <- read_csv(file.path(A, "simulation_strong_ALLGENE_PROBE_SUMMARY.csv"))
combined <- rbind(base_probe, strong_probe)
checks[[length(checks) + 1L]] <- rec("exactly_three_predeclared_allgene_examples", nrow(base_probe) == 2L && nrow(strong_probe) == 1L &&
  setequal(base_probe$scenario, c("overlap", "separation")) && identical(strong_probe$scenario, "strong") && all(combined$replicate == 1L),
  sprintf("existing=%d; strong=%d; total=%d; statuses=%s", nrow(base_probe), nrow(strong_probe), nrow(combined), paste(combined$status, collapse = "/")))

result <- do.call(rbind, checks)
write.csv(result, file.path(A, "simulation_strong_VALIDATION.csv"), row.names = FALSE, na = "", quote = TRUE)
resources <- read_csv(file.path(A, "simulation_strong_repeat_resource_metrics.csv"))
cat("Independent strong checks=", nrow(result), "; passed=", sum(result$pass), "; failed=", sum(!result$pass),
    "; estimated VIF min/median/max=", sprintf("%.3f/%.3f/%.3f", min(vif), median(vif), max(vif)),
    "; gate failures=", gate_fail, "; observed R GC/vector heap max across repeats MB=", sprintf("%.1f", max(resources$gc_peak_heap_mb_observed)),
    "; warning/failure records=", nrow(warnings), ".\n", sep = "")
if (any(!result$pass)) quit(save = "no", status = 1L)
