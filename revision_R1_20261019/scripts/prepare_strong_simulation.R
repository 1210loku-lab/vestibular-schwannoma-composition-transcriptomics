#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
ROOT <- Sys.getenv("VS_INPUT_ROOT", unset = ".")
EXPORT_ROOT <- Sys.getenv("VS_OUTPUT_ROOT", unset = ".")
A <- file.path(EXPORT_ROOT, "outputs")
SCRIPT <- file.path(EXPORT_ROOT, "scripts", "prepare_strong_simulation.R")
RUNNER_SCRIPT <- file.path(EXPORT_ROOT, "scripts", "run_strong_simulation.R")
PROBE_PLAN <- file.path(A, "simulation_strong_ALLGENE_PROBE_PLAN.md")
PROBE_SCRIPT <- file.path(EXPORT_ROOT, "scripts", "probe_all_gene_strong_simulation.R")
VALIDATOR_SCRIPT <- file.path(EXPORT_ROOT, "scripts", "validate_strong_simulation.R")
cell_types <- c("Schwann/肿瘤细胞", "巨噬/髓系", "T细胞", "成纤维", "内皮", "周细胞")
features <- c("Schwann", "Myeloid", "T_cell", "Fibroblast", "Endothelial", "Pericyte")
N_REPEAT <- 20L
N_CONTROL <- 9L
N_TUMOR <- 31L
CONCENTRATION <- 100
SEEDS <- 2026101801:2026101820
NOISE_SD_LOG2 <- 0.12
sha256 <- function(path) {
  z <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  if (length(z) != 1L || !grepl("^[0-9a-fA-F]{64}\\s", z)) stop("SHA-256 failed for ", path)
  strsplit(z, "\\s+")[[1]][1]
}
read_csv <- function(path) utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
read_gz <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  utils::read.csv(con, check.names = FALSE, stringsAsFactors = FALSE)
}
write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE, na = "", quote = TRUE)
write_gz <- function(x, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.csv(x, con, row.names = FALSE, na = "", quote = TRUE)
}

message("Strong extension preparation started; no estimator is called.")
strong_md <- file.path(A, "simulation_strong_PROTOCOL.md")
strong_json <- file.path(A, "simulation_strong_PROTOCOL.json")
base_manifest <- file.path(A, "simulation_FROZEN_INPUT_MANIFEST.csv")
if (!file.exists(strong_md) || !file.exists(strong_json) || !file.exists(base_manifest)) stop("A frozen protocol or base manifest is missing")
protocol <- jsonlite::fromJSON(strong_json, simplifyVector = FALSE)
if (!identical(as.numeric(unlist(protocol$scope$new_seed_range_inclusive)), as.numeric(c(2026101801, 2026101820))) ||
    protocol$scope$new_concentration != CONCENTRATION || protocol$scope$new_repeats != N_REPEAT) {
  stop("Strong extension parameters differ from the frozen protocol")
}
frozen_time <- protocol$frozen_before_estimator_outcomes
base <- read_csv(base_manifest)
base_hashes <- vapply(base$path, sha256, character(1))
if (any(base_hashes != base$sha256)) stop("A corrected base simulation frozen-input hash changed")
prior_probe_summary_path <- file.path(A, "simulation_ALLGENE_PROBE_SUMMARY.csv")
prior_probe_manifest_path <- file.path(A, "simulation_ALLGENE_PROBE_OUTPUT_MANIFEST.csv")
prior_probe_input_manifest_path <- file.path(A, "simulation_ALLGENE_PROBE_INPUT_MANIFEST.csv")
prior_probe <- read_csv(prior_probe_summary_path)
prior_probe_manifest <- read_csv(prior_probe_manifest_path)
prior_probe_hashes <- vapply(prior_probe_manifest$path, sha256, character(1))
if (any(prior_probe_hashes != prior_probe_manifest$sha256) || nrow(prior_probe) != 2L ||
    !setequal(prior_probe$scenario, c("overlap", "separation")) ||
    any(prior_probe$replicate != 1L) || any(prior_probe$status != "ok")) {
  stop("The two predeclared overlap/moderate all-gene examples failed their frozen integrity/status check")
}

truth <- read_csv(file.path(A, "simulation_gene_truth_and_marker_allocation.csv"))
profile_tab <- read_gz(file.path(A, "simulation_baseline_linear_profiles.csv.gz"))
profile <- as.matrix(profile_tab[, -1, drop = FALSE])
storage.mode(profile) <- "double"
rownames(profile) <- profile_tab[[1]]
if (nrow(truth) != 1500L || !identical(rownames(profile), truth$gene) || ncol(profile) != 6L ||
    sum(truth$reserved_marker) != 300L) stop("Reused truth/profile inputs do not match frozen simulation design")
beta <- setNames(truth$beta, truth$gene)
if (any(!is.finite(beta)) || any(abs(beta) > 1.5)) stop("Unexpected beta in shared truth")
old_meta <- read_csv(file.path(A, "simulation_mixture_sample_metadata.csv"))
if (length(intersect(unique(as.integer(old_meta$seed)), SEEDS))) stop("Strong seeds overlap earlier simulation mixtures")

bulk_path <- file.path(A, "bulk_composition_distribution.csv")
bulk_composition <- read_csv(bulk_path)
bulk_composition <- bulk_composition[bulk_composition$reference == "MuSiC" & bulk_composition$feature %in% features & bulk_composition$group %in% c("Control", "Tumor"), , drop = FALSE]
get_bulk_mean <- function(group) {
  d <- bulk_composition[bulk_composition$group == group, , drop = FALSE]
  vals <- d$mean[match(features, d$feature)]
  if (length(vals) != length(features) || anyNA(vals) || any(!is.finite(vals)) || any(vals < 0) || sum(vals) <= 0) {
    stop("Invalid bulk-model group means: ", group)
  }
  vals
}
bulk_control_mean <- get_bulk_mean("Control")
bulk_tumor_mean <- get_bulk_mean("Tumor")
control_center <- bulk_control_mean / sum(bulk_control_mean)
tumor_center <- bulk_tumor_mean / sum(bulk_tumor_mean)
base_centers <- read_csv(file.path(A, "simulation_scenario_centers_and_bulk_contrast.csv"))
moderate <- base_centers[base_centers$scenario == "separation", , drop = FALSE]
if (nrow(moderate) != length(features) ||
    max(abs(moderate$bulk_control_normalized - control_center)) > 1e-12 ||
    max(abs(moderate$bulk_tumor_normalized - tumor_center)) > 1e-12 ||
    max(abs(moderate$control_center - control_center)) > 1e-12 ||
    max(abs(moderate$tumor_center - tumor_center)) > 1e-12) {
  stop("Strong scenario centers do not exactly reuse the corrected bulk-model-derived centers")
}
if (max(abs(control_center - tumor_center)) <= 0.10) stop("bulk-model Control/Tumor centers are not distinctly separated")

mixture_y <- matrix(NA_real_, nrow = nrow(profile), ncol = N_REPEAT * (N_CONTROL + N_TUMOR),
                    dimnames = list(rownames(profile), NULL))
sample_rows <- vector("list", ncol(mixture_y))
truth_rows <- vector("list", ncol(mixture_y) * length(features))
sample_i <- 0L
truth_i <- 0L
for (rep_i in seq_len(N_REPEAT)) {
  seed <- SEEDS[rep_i]
  set.seed(seed)
  group_vec <- c(rep("Control", N_CONTROL), rep("Tumor", N_TUMOR))
  for (within_i in seq_along(group_vec)) {
    group <- group_vec[within_i]
    center <- if (group == "Control") control_center else tumor_center
    p <- stats::rgamma(length(center), shape = center * CONCENTRATION, rate = 1)
    if (any(!is.finite(p)) || sum(p) <= 0) stop("Invalid strong Dirichlet composition draw")
    p <- p / sum(p)
    sid <- sprintf("STR_R%02d_%s_%02d", rep_i, group, if (group == "Control") within_i else within_i - N_CONTROL)
    tumor <- as.integer(group == "Tumor")
    fold <- 2^(unname(beta[rownames(profile)]) * tumor)
    linear <- drop(profile %*% p) * fold
    if (any(!is.finite(linear)) || any(linear <= 0)) stop("Invalid strong linear mixture")
    y <- log2(linear) + stats::rnorm(length(linear), 0, NOISE_SD_LOG2)
    sample_i <- sample_i + 1L
    mixture_y[, sample_i] <- y
    sample_rows[[sample_i]] <- data.frame(
      scenario = "strong", concentration = CONCENTRATION, replicate = rep_i,
      seed = seed, sample = sid, group = group, dirichlet_center_group = group,
      stringsAsFactors = FALSE
    )
    rna_mass <- p * colSums(sweep(profile, 1, fold, `*`))
    q <- rna_mass / sum(rna_mass)
    for (ct_i in seq_along(features)) {
      truth_i <- truth_i + 1L
      truth_rows[[truth_i]] <- data.frame(
        scenario = "strong", concentration = CONCENTRATION, replicate = rep_i,
        seed = seed, sample = sid, group = group, dirichlet_center_group = group,
        feature = features[ct_i], cell_type = cell_types[ct_i],
        assigned_center_component = center[ct_i], cell_fraction = p[ct_i],
        rna_contribution_fraction = q[ct_i], stringsAsFactors = FALSE
      )
    }
  }
}
sample_meta <- do.call(rbind, sample_rows)
fraction_truth <- do.call(rbind, truth_rows)
colnames(mixture_y) <- sample_meta$sample
if (nrow(sample_meta) != 800L || nrow(fraction_truth) != 4800L || ncol(mixture_y) != 800L ||
    length(unique(sample_meta$seed)) != 20L || !identical(sort(unique(sample_meta$seed)), SEEDS)) {
  stop("Strong extension sample or seed count differs from fixed scope")
}
if (length(intersect(unique(as.integer(old_meta$seed)), unique(sample_meta$seed))) != 0L) stop("Strong seeds overlap prior samples")
if (any(sample_meta$group != sample_meta$dirichlet_center_group) ||
    any(fraction_truth$group != fraction_truth$dirichlet_center_group)) stop("A sample was assigned the wrong group center")
sum_cell <- aggregate(cell_fraction ~ sample, fraction_truth, sum)
sum_rna <- aggregate(rna_contribution_fraction ~ sample, fraction_truth, sum)
if (max(abs(sum_cell$cell_fraction - 1)) > 1e-10 || max(abs(sum_rna$rna_contribution_fraction - 1)) > 1e-10) {
  stop("Strong cell/RNA truth fractions do not sum to one")
}

group_means <- aggregate(cbind(cell_fraction, rna_contribution_fraction) ~ group + feature,
                         fraction_truth, mean)
center_truth <- data.frame(
  scenario = "strong", concentration = CONCENTRATION,
  feature = features, cell_type = cell_types,
  bulk_control_mean = bulk_control_mean, bulk_tumor_mean = bulk_tumor_mean,
  control_center = control_center, tumor_center = tumor_center,
  center_tumor_minus_control = tumor_center - control_center,
  stringsAsFactors = FALSE
)
for (group in c("Control", "Tumor")) {
  gm <- group_means[group_means$group == group, , drop = FALSE]
  gm <- gm[match(features, gm$feature), , drop = FALSE]
  center_truth[[paste0(tolower(group), "_realized_mean_cell_fraction")]] <- gm$cell_fraction
  center_truth[[paste0(tolower(group), "_realized_mean_rna_contribution")]] <- gm$rna_contribution_fraction
}
center_truth$realized_mean_cell_fraction_tumor_minus_control <-
  center_truth$tumor_realized_mean_cell_fraction - center_truth$control_realized_mean_cell_fraction
center_truth$realized_mean_rna_contribution_tumor_minus_control <-
  center_truth$tumor_realized_mean_rna_contribution - center_truth$control_realized_mean_rna_contribution
if (max(abs(center_truth$realized_mean_cell_fraction_tumor_minus_control)) < 0.10) {
  stop("Strong realized group composition contrast did not pass the pre-fit design check")
}

# Load-bearing same-composition, zero-noise identity across every gene.
p_test <- as.numeric(fraction_truth[fraction_truth$sample == sample_meta$sample[1], "cell_fraction"])
linear_control <- drop(profile %*% p_test)
linear_tumor <- drop(profile %*% p_test) * 2^unname(beta[rownames(profile)])
identity <- data.frame(
  gene = rownames(profile), beta = unname(beta[rownames(profile)]),
  noiseless_fixed_composition_log2_contrast = log2(linear_tumor / linear_control),
  absolute_error = abs(log2(linear_tumor / linear_control) - unname(beta[rownames(profile)])),
  stringsAsFactors = FALSE
)
identity_max <- max(identity$absolute_error)
if (identity_max > 1e-12) stop("Strong fixed-composition noiseless contrast does not equal beta")

write_gz(data.frame(gene = rownames(mixture_y), mixture_y, check.names = FALSE),
         file.path(A, "simulation_strong_mixture_log2_expression.csv.gz"))
write_csv(sample_meta, file.path(A, "simulation_strong_sample_metadata.csv"))
write_gz(fraction_truth, file.path(A, "simulation_strong_fraction_truth.csv.gz"))
write_csv(identity, file.path(A, "simulation_strong_noiseless_truth_identity_check.csv"))
write_csv(center_truth, file.path(A, "simulation_strong_centers_and_truth_contrast.csv"))

checks <- data.frame(
  check = c("new_seed_count_exact", "strong_seed_nonoverlap", "group_center_alignment_all_samples",
            "centers_exactly_match_bulk_and_moderate", "bulk_centers_are_distinct", "cell_and_RNA_truth_sum_to_one",
            "realized_group_truth_separation", "noiseless_fixed_composition_beta_identity"),
  pass = c(length(unique(sample_meta$seed)) == 20L,
           length(intersect(unique(as.integer(old_meta$seed)), unique(sample_meta$seed))) == 0L,
           all(sample_meta$group == sample_meta$dirichlet_center_group) &&
             all(fraction_truth$group == fraction_truth$dirichlet_center_group),
           max(abs(center_truth$control_center - moderate$control_center)) <= 1e-12 &&
             max(abs(center_truth$tumor_center - moderate$tumor_center)) <= 1e-12,
           max(abs(center_truth$center_tumor_minus_control)) > 0.10,
           max(abs(sum_cell$cell_fraction - 1)) < 1e-10 && max(abs(sum_rna$rna_contribution_fraction - 1)) < 1e-10,
           max(abs(center_truth$realized_mean_cell_fraction_tumor_minus_control)) > 0.10,
           identity_max < 1e-12),
  detail = c(sprintf("%d unique seeds; range %d-%d", length(unique(sample_meta$seed)), min(sample_meta$seed), max(sample_meta$seed)),
             sprintf("overlap with old mixture seeds=%d", length(intersect(unique(as.integer(old_meta$seed)), unique(sample_meta$seed)))),
             sprintf("samples=%d; rows aligned to their own group", nrow(sample_meta)),
             sprintf("max center discrepancy=%.3g", max(c(abs(center_truth$control_center - moderate$control_center),
                                                            abs(center_truth$tumor_center - moderate$tumor_center)))),
             sprintf("max component contrast=%.6f", max(abs(center_truth$center_tumor_minus_control))),
             sprintf("cell sum error=%.3g; RNA sum error=%.3g", max(abs(sum_cell$cell_fraction - 1)), max(abs(sum_rna$rna_contribution_fraction - 1))),
             sprintf("max realized mean cell-fraction contrast=%.6f; RNA=%.6f",
                     max(abs(center_truth$realized_mean_cell_fraction_tumor_minus_control)),
                     max(abs(center_truth$realized_mean_rna_contribution_tumor_minus_control))),
             sprintf("max absolute error across %d genes=%.3g", nrow(identity), identity_max)),
  stringsAsFactors = FALSE
)
if (any(!checks$pass)) stop("A strong pre-fit truth/design gate failed")
write_csv(checks, file.path(A, "simulation_strong_DESIGN_AND_TRUTH_GATE.csv"))

frozen_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
prior_probe_paths <- unique(c(prior_probe_manifest$path, prior_probe_manifest_path,
                              prior_probe_input_manifest_path, prior_probe_summary_path))
input_paths <- c(
  strong_md, strong_json, SCRIPT,
  RUNNER_SCRIPT, PROBE_PLAN, PROBE_SCRIPT, VALIDATOR_SCRIPT,
  file.path(A, "simulation_PILOT_PROTOCOL.md"), file.path(A, "simulation_PILOT_PROTOCOL.json"),
  base_manifest, file.path(A, "simulation_gene_truth_and_marker_allocation.csv"),
  file.path(A, "simulation_baseline_linear_profiles.csv.gz"),
  file.path(A, "simulation_synthetic_reference_counts.csv.gz"),
  file.path(A, "simulation_synthetic_reference_cell_metadata.csv"),
  file.path(A, "bulk_composition_distribution.csv"),
  file.path(A, "simulation_scenario_centers_and_bulk_contrast.csv"),
  file.path(A, "simulation_mixture_sample_metadata.csv"),
  file.path(A, "simulation_strong_mixture_log2_expression.csv.gz"),
  file.path(A, "simulation_strong_sample_metadata.csv"),
  file.path(A, "simulation_strong_fraction_truth.csv.gz"),
  file.path(A, "simulation_strong_noiseless_truth_identity_check.csv"),
  file.path(A, "simulation_strong_centers_and_truth_contrast.csv"),
  file.path(A, "simulation_strong_DESIGN_AND_TRUTH_GATE.csv"),
  prior_probe_paths
)
manifest <- data.frame(
  path = normalizePath(input_paths, winslash = "/", mustWork = TRUE),
  role = c("strong_protocol_md", "strong_protocol_json", "strong_generator_script",
           "strong_primary_runner_script", "strong_allgene_probe_plan", "strong_allgene_probe_script",
           "independent_validator_script",
           "base_protocol_md", "base_protocol_json", "base_corrected_manifest",
           "reused_gene_truth", "reused_linear_profile", "reused_synthetic_reference_counts",
           "reused_synthetic_reference_metadata", "bulk_center_source", "reused_corrected_bulk_centers",
           "old_seed_overlap_check", "new_strong_mixture", "new_strong_sample_metadata",
           "new_strong_cell_and_RNA_truth", "new_strong_identity_check", "new_strong_center_truth",
           "pre_fit_design_truth_gate", rep("preexisting_overlap_and_moderate_allgene_artifact", length(prior_probe_paths))),
  bytes = file.info(input_paths)$size,
  sha256 = vapply(input_paths, sha256, character(1)),
  freeze_timestamp_local = frozen_at,
  stringsAsFactors = FALSE
)
write_csv(manifest, file.path(A, "simulation_strong_FROZEN_INPUT_MANIFEST.csv"))
cat("Strong extension frozen before estimator outcomes at ", frozen_at, ".\n", sep = "")
cat("protocol SHA256 md/json=", sha256(strong_md), "/", sha256(strong_json), "\n", sep = "")
cat("unique new seeds=", length(unique(sample_meta$seed)), "; samples=", nrow(sample_meta),
    "; truth rows=", nrow(fraction_truth), "; exact bulk-model center max difference=",
    format(max(abs(center_truth$center_tumor_minus_control)), digits = 7),
    "; realized cell contrast max=",
    format(max(abs(center_truth$realized_mean_cell_fraction_tumor_minus_control)), digits = 7),
    "; identity max error=", format(identity_max, scientific = TRUE), "\n", sep = "")
cat("pre-fit checks passed=", sum(checks$pass), "/", nrow(checks), "; frozen input records=", nrow(manifest), "\n", sep = "")
