#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

ROOT <- Sys.getenv("VS_INPUT_ROOT", unset = ".")
EXPORT_ROOT <- Sys.getenv("VS_OUTPUT_ROOT", unset = ".")
ANALYSIS <- file.path(EXPORT_ROOT, "outputs")
MEAN_PATH <- file.path(ROOT, "data", "processed", "GSE230375_celltype_meanexpr.csv")
bulk_COMP_PATH <- file.path(ANALYSIS, "bulk_composition_distribution.csv")
PROTOCOL_MD <- file.path(ANALYSIS, "simulation_PILOT_PROTOCOL.md")
PROTOCOL_JSON <- file.path(ANALYSIS, "simulation_PILOT_PROTOCOL.json")
SCRIPT_PATH <- file.path(EXPORT_ROOT, "scripts", "prepare_controlled_simulation.R")

cell_types <- c("Schwann/肿瘤细胞", "巨噬/髓系", "T细胞", "成纤维", "内皮", "周细胞")
features <- c("Schwann", "Myeloid", "T_cell", "Fibroblast", "Endothelial", "Pericyte")
N_GENES <- 1500L
N_MARKERS_PER_CT <- 50L
N_TARGET_PER_CLASS <- 100L
N_BACKGROUND <- 800L
N_DONORS <- 5L
N_CELLS_PER_DONOR_CT <- 30L
BASELINE_COUNTS_PER_CELL <- 10000
DONOR_SD_LOG <- 0.10
NB_SIZE <- 20
N_REPEATS <- 20L
N_CONTROL <- 9L
N_TUMOR <- 31L
NOISE_SD_LOG2 <- 0.12
ALLOCATION_SEED <- 20261015L
REFERENCE_SEED <- 2026101510L

write_csv <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, na = "", quote = TRUE)
}

write_gz_csv <- function(x, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.csv(x, con, row.names = FALSE, na = "", quote = TRUE)
}

sha256 <- function(path) {
  out <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  if (length(out) != 1L || !grepl("^[0-9a-fA-F]{64}\\s", out)) {
    stop("Could not compute SHA-256 for ", path, ": ", paste(out, collapse = " | "))
  }
  strsplit(out, "\\s+")[[1]][1]
}

message("simulation input preparation started; no deconvolution or outcome fit is run in this script.")
if (!file.exists(PROTOCOL_MD) || !file.exists(PROTOCOL_JSON)) stop("Frozen protocol files are missing")
if (!file.exists(MEAN_PATH) || !file.exists(bulk_COMP_PATH)) stop("A declared input is missing")
protocol_txt <- paste(readLines(PROTOCOL_MD, warn = FALSE), collapse = "\n")
if (!grepl("Frozen before estimator outcomes", protocol_txt, fixed = TRUE)) stop("Protocol is not marked frozen")

mean_tab <- utils::read.csv(MEAN_PATH, check.names = FALSE)
if (!identical(names(mean_tab)[-1], cell_types)) stop("Unexpected cell type order in mean-expression input")
genes_all <- as.character(mean_tab[[1]])
if (anyDuplicated(genes_all) || anyNA(genes_all) || any(!nzchar(genes_all))) stop("Invalid gene IDs in background")
mean_log <- as.matrix(mean_tab[, -1, drop = FALSE])
storage.mode(mean_log) <- "double"
rownames(mean_log) <- genes_all
if (any(!is.finite(mean_log)) || any(mean_log < 0)) stop("Background mean expression must be finite and nonnegative")
source_linear <- expm1(mean_log)
rownames(source_linear) <- genes_all
eligible <- rownames(source_linear)[rowSums(source_linear) >= 2 & apply(source_linear, 1, max) >= 0.5]
if (length(eligible) < N_GENES) stop("Too few expressed genes for fixed allocation")

# Freeze the marker set first. Targets are selected only from genes outside this set.
max_other <- sapply(seq_along(cell_types), function(j) {
  apply(source_linear[, -j, drop = FALSE], 1, max)
})
specificity <- log2((source_linear + 0.05) / (max_other + 0.05))
rownames(specificity) <- genes_all
markers <- character()
marker_assignment <- character()
for (j in seq_along(cell_types)) {
  candidates <- eligible[order(-specificity[eligible, j], -source_linear[eligible, j], eligible)]
  candidates <- candidates[!candidates %in% markers]
  picked <- head(candidates, N_MARKERS_PER_CT)
  if (length(picked) != N_MARKERS_PER_CT) stop("Could not allocate 50 unique markers for ", cell_types[j])
  markers <- c(markers, picked)
  marker_assignment <- c(marker_assignment, rep(cell_types[j], length(picked)))
}

set.seed(ALLOCATION_SEED)
remaining <- setdiff(eligible, markers)
heterogeneity <- apply(log2(source_linear[remaining, , drop = FALSE] + 0.05), 1, stats::sd)
heterogeneous_pool <- remaining[order(-heterogeneity[remaining], remaining)]
heterogeneous_targets <- head(heterogeneous_pool, 2L * N_TARGET_PER_CLASS)
if (length(heterogeneous_targets) != 2L * N_TARGET_PER_CLASS) stop("Insufficient heterogeneous target genes")
composition_only <- heterogeneous_targets[seq_len(N_TARGET_PER_CLASS)]
both <- heterogeneous_targets[N_TARGET_PER_CLASS + seq_len(N_TARGET_PER_CLASS)]
remaining <- setdiff(remaining, heterogeneous_targets)
homogeneous_pool <- remaining[order(-rowMeans(source_linear[remaining, , drop = FALSE]), remaining)]
homogeneous_pool <- head(homogeneous_pool, 1000L)
homogeneous_targets <- sample(homogeneous_pool, 2L * N_TARGET_PER_CLASS, replace = FALSE)
null_target <- homogeneous_targets[seq_len(N_TARGET_PER_CLASS)]
expression_only <- homogeneous_targets[N_TARGET_PER_CLASS + seq_len(N_TARGET_PER_CLASS)]
remaining <- setdiff(remaining, homogeneous_targets)
background <- sample(remaining, N_BACKGROUND, replace = FALSE)
gene_ids <- c(markers, null_target, composition_only, expression_only, both, background)
if (length(gene_ids) != N_GENES || anyDuplicated(gene_ids)) stop("Frozen gene allocation is not exactly 1,500 unique genes")

truth_class <- setNames(rep("non_target_background", length(gene_ids)), gene_ids)
truth_class[markers] <- "reserved_marker"
truth_class[null_target] <- "null"
truth_class[composition_only] <- "composition_only"
truth_class[expression_only] <- "expression_only"
truth_class[both] <- "both"
beta <- setNames(rep(0, length(gene_ids)), gene_ids)
for (class_genes in list(expression_only, both)) {
  effect_values <- sample(rep(c(-1.5, 1.5), each = 50L), 100L, replace = FALSE)
  beta[class_genes] <- effect_values
}

# Build positive linear profiles. Homogeneous target profiles are set equal after
# reference-source scaling; remaining genes are rescaled to preserve 10,000 total.
profile <- source_linear[gene_ids, cell_types, drop = FALSE]
profile[!is.finite(profile) | profile < 0] <- 0
profile <- pmax(profile, 0.01)
homogeneous_genes <- c(null_target, expression_only)
homogeneous_value <- pmax(rowMeans(source_linear[homogeneous_genes, cell_types, drop = FALSE]), 0.1) * 0.25
names(homogeneous_value) <- homogeneous_genes
profile[homogeneous_genes, ] <- homogeneous_value[homogeneous_genes]
nonhomogeneous <- setdiff(gene_ids, homogeneous_genes)
homogeneous_total <- sum(homogeneous_value)
if (homogeneous_total >= BASELINE_COUNTS_PER_CELL) stop("Homogeneous target totals exceed the fixed cell total")
for (ct in cell_types) {
  denom <- sum(profile[nonhomogeneous, ct])
  profile[nonhomogeneous, ct] <- profile[nonhomogeneous, ct] *
    ((BASELINE_COUNTS_PER_CELL - homogeneous_total) / denom)
  profile[homogeneous_genes, ct] <- homogeneous_value[homogeneous_genes]
}
rownames(profile) <- gene_ids
if (any(!is.finite(profile)) || any(profile <= 0)) stop("Baseline profiles must be positive and finite")
if (max(abs(colSums(profile) - BASELINE_COUNTS_PER_CELL)) > 1e-8) stop("Baseline profile totals are not equal")
if (max(abs(apply(profile[null_target, , drop = FALSE], 1, stats::sd))) > 1e-12 ||
    max(abs(apply(profile[expression_only, , drop = FALSE], 1, stats::sd))) > 1e-12) {
  stop("Null and expression-only target profiles must be constant across compartments")
}

gene_truth <- data.frame(
  gene = gene_ids,
  gene_class = unname(truth_class[gene_ids]),
  beta = unname(beta[gene_ids]),
  reserved_marker = gene_ids %in% markers,
  marker_compartment = c(marker_assignment, rep(NA_character_, N_GENES - length(markers))),
  source_mean_linear = rowMeans(source_linear[gene_ids, cell_types, drop = FALSE]),
  background_heterogeneity_sd_log2 = apply(log2(source_linear[gene_ids, cell_types, drop = FALSE] + 0.05), 1, stats::sd),
  stringsAsFactors = FALSE
)
write_csv(gene_truth, file.path(ANALYSIS, "simulation_gene_truth_and_marker_allocation.csv"))
write_gz_csv(data.frame(gene = gene_ids, profile, check.names = FALSE),
             file.path(ANALYSIS, "simulation_baseline_linear_profiles.csv.gz"))

# Synthetic reference donors/cells are generated independently of mixture draws.
set.seed(REFERENCE_SEED)
reference_counts <- matrix(0L, nrow = N_GENES,
                           ncol = N_DONORS * length(cell_types) * N_CELLS_PER_DONOR_CT,
                           dimnames = list(gene_ids, NULL))
reference_meta <- vector("list", ncol(reference_counts))
cell_col <- 0L
for (d in seq_len(N_DONORS)) {
  donor <- sprintf("SIM_D%02d", d)
  for (ct in cell_types) {
    donor_profile <- profile[, ct] * exp(stats::rnorm(N_GENES, 0, DONOR_SD_LOG))
    donor_profile <- donor_profile * (BASELINE_COUNTS_PER_CELL / sum(donor_profile))
    for (cell_i in seq_len(N_CELLS_PER_DONOR_CT)) {
      cell_col <- cell_col + 1L
      reference_counts[, cell_col] <- as.integer(stats::rnbinom(N_GENES, size = NB_SIZE, mu = donor_profile))
      reference_meta[[cell_col]] <- data.frame(
        cell = sprintf("%s_CT%02d_C%02d", donor, match(ct, cell_types), cell_i),
        sample = donor,
        cell_type = ct,
        stringsAsFactors = FALSE
      )
    }
  }
}
colnames(reference_counts) <- vapply(reference_meta, `[[`, character(1), "cell")
reference_meta <- do.call(rbind, reference_meta)
write_gz_csv(data.frame(gene = gene_ids, reference_counts, check.names = FALSE),
             file.path(ANALYSIS, "simulation_synthetic_reference_counts.csv.gz"))
write_csv(reference_meta, file.path(ANALYSIS, "simulation_synthetic_reference_cell_metadata.csv"))

# Group centers for the separated scenario are imported as design parameters only.
bulk_composition <- utils::read.csv(bulk_COMP_PATH, stringsAsFactors = FALSE)
bulk_music <- bulk_composition[bulk_composition$reference == "MuSiC" & bulk_composition$feature %in% features & bulk_composition$group %in% c("Control", "Tumor"), ]
if (nrow(bulk_music) != 12L) stop("bulk-model MuSiC composition centers do not contain 12 expected rows")
values_from_bulk <- function(group_name) {
  d <- bulk_music[bulk_music$group == group_name, , drop = FALSE]
  out <- d$mean[match(features, d$feature)]
  if (anyNA(out) || any(!is.finite(out)) || any(out < 0) || sum(out) <= 0) {
    stop("Missing or invalid bulk-model feature center for ", group_name)
  }
  out
}
bulk_control_mean <- values_from_bulk("Control")
bulk_tumor_mean <- values_from_bulk("Tumor")
center_from_bulk <- function(group_name) {
  out <- values_from_bulk(group_name)
  out / sum(out)
}
bulk_control_center <- center_from_bulk("Control")
bulk_tumor_center <- center_from_bulk("Tumor")
if (anyNA(bulk_control_mean) || anyNA(bulk_tumor_mean) ||
    max(abs(bulk_control_center - bulk_tumor_center)) < 0.10) {
  stop("bulk-model Control/Tumor composition centers are missing or insufficiently distinct")
}
scenario_defs <- list(
  overlap = list(
    scenario_index = 1L,
    control_center = c(0.32, 0.15, 0.12, 0.14, 0.17, 0.10),
    tumor_center = c(0.36, 0.15, 0.12, 0.13, 0.16, 0.08),
    concentration = 50
  ),
  separation = list(
    scenario_index = 2L,
    control_center = center_from_bulk("Control"),
    tumor_center = center_from_bulk("Tumor"),
    concentration = 50
  )
)
if (max(abs(scenario_defs$separation$control_center - bulk_control_center)) > 1e-12 ||
    max(abs(scenario_defs$separation$tumor_center - bulk_tumor_center)) > 1e-12) {
  stop("Separated-scenario centers do not exactly match normalized bulk-model group means")
}
if (any(abs(vapply(scenario_defs, function(z) sum(z$control_center), numeric(1)) - 1) > 1e-10) ||
    any(abs(vapply(scenario_defs, function(z) sum(z$tumor_center), numeric(1)) - 1) > 1e-10)) {
  stop("Scenario composition centers must sum to one")
}
draw_dirichlet <- function(center, concentration) {
  value <- stats::rgamma(length(center), shape = center * concentration, rate = 1)
  if (any(!is.finite(value)) || sum(value) <= 0) stop("Invalid Dirichlet draw")
  value / sum(value)
}

n_total <- length(scenario_defs) * N_REPEATS * (N_CONTROL + N_TUMOR)
mixture_y <- matrix(NA_real_, nrow = N_GENES, ncol = n_total,
                    dimnames = list(gene_ids, NULL))
mixture_sample_ids <- character(n_total)
sample_meta <- vector("list", n_total)
composition_truth <- vector("list", n_total * length(cell_types))
meta_i <- 0L
comp_i <- 0L
for (scenario_name in names(scenario_defs)) {
  def <- scenario_defs[[scenario_name]]
  for (rep_i in seq_len(N_REPEATS)) {
    seed <- 2026101500L + 100L * def$scenario_index + rep_i
    set.seed(seed)
    group_vec <- c(rep("Control", N_CONTROL), rep("Tumor", N_TUMOR))
    for (sample_i in seq_along(group_vec)) {
      group <- group_vec[sample_i]
      center <- if (group == "Control") def$control_center else def$tumor_center
      p <- draw_dirichlet(center, def$concentration)
      sample_id <- sprintf("%s_R%02d_%s_%02d", toupper(substr(scenario_name, 1, 2)), rep_i,
                           group, if (group == "Control") sample_i else sample_i - N_CONTROL)
      tumor <- as.integer(group == "Tumor")
      fold <- 2^(unname(beta) * tumor)
      linear <- drop(profile %*% p) * fold
      if (any(!is.finite(linear)) || any(linear <= 0)) stop("Invalid generated linear mixture")
      y <- log2(linear) + stats::rnorm(N_GENES, 0, NOISE_SD_LOG2)
      meta_i <- meta_i + 1L
      mixture_y[, meta_i] <- y
      mixture_sample_ids[meta_i] <- sample_id
      sample_meta[[meta_i]] <- data.frame(
        scenario = scenario_name,
        scenario_index = def$scenario_index,
        replicate = rep_i,
        seed = seed,
        sample = sample_id,
        group = group,
        stringsAsFactors = FALSE
      )
      rna_mass <- p * colSums(sweep(profile, 1, fold, `*`))
      rna_contribution <- rna_mass / sum(rna_mass)
      for (ct_i in seq_along(cell_types)) {
        comp_i <- comp_i + 1L
        composition_truth[[comp_i]] <- data.frame(
          scenario = scenario_name,
          replicate = rep_i,
          sample = sample_id,
          group = group,
          cell_type = cell_types[ct_i],
          cell_fraction = p[ct_i],
          rna_contribution_fraction = rna_contribution[ct_i],
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
colnames(mixture_y) <- mixture_sample_ids
sample_meta <- do.call(rbind, sample_meta)
composition_truth <- do.call(rbind, composition_truth)
if (ncol(mixture_y) != n_total || nrow(sample_meta) != n_total || nrow(composition_truth) != n_total * 6L) {
  stop("Unexpected frozen mixture count")
}
truth_group_means <- aggregate(
  cbind(cell_fraction, rna_contribution_fraction) ~ scenario + group + cell_type,
  data = composition_truth, FUN = mean
)
center_comparison <- do.call(rbind, lapply(names(scenario_defs), function(scenario_name) {
  def <- scenario_defs[[scenario_name]]
  out <- data.frame(
    scenario = scenario_name,
    feature = features,
    cell_type = cell_types,
    control_center = def$control_center,
    tumor_center = def$tumor_center,
    center_tumor_minus_control = def$tumor_center - def$control_center,
    bulk_control_mean = if (scenario_name == "separation") bulk_control_mean else NA_real_,
    bulk_tumor_mean = if (scenario_name == "separation") bulk_tumor_mean else NA_real_,
    bulk_control_normalized = if (scenario_name == "separation") bulk_control_center else NA_real_,
    bulk_tumor_normalized = if (scenario_name == "separation") bulk_tumor_center else NA_real_,
    stringsAsFactors = FALSE
  )
  for (group_name in c("Control", "Tumor")) {
    rows <- truth_group_means[truth_group_means$scenario == scenario_name &
                                truth_group_means$group == group_name, , drop = FALSE]
    rows <- rows[match(cell_types, rows$cell_type), , drop = FALSE]
    out[[paste0(tolower(group_name), "_realized_mean_cell_fraction")]] <- rows$cell_fraction
    out[[paste0(tolower(group_name), "_realized_mean_rna_contribution")]] <- rows$rna_contribution_fraction
  }
  out$realized_mean_cell_fraction_tumor_minus_control <-
    out$tumor_realized_mean_cell_fraction - out$control_realized_mean_cell_fraction
  out$realized_mean_rna_contribution_tumor_minus_control <-
    out$tumor_realized_mean_rna_contribution - out$control_realized_mean_rna_contribution
  out
}))
sep_truth_delta <- center_comparison$realized_mean_cell_fraction_tumor_minus_control[
  center_comparison$scenario == "separation"
]
if (max(abs(sep_truth_delta)) < 0.10 ||
    max(abs(center_comparison$center_tumor_minus_control[
      center_comparison$scenario == "separation"
    ])) < 0.10) {
  stop("Separated scenario truth does not retain a distinct Control/Tumor composition contrast")
}
write_csv(center_comparison, file.path(ANALYSIS, "simulation_scenario_centers_and_bulk_contrast.csv"))
if (any(abs(aggregate(cell_fraction ~ scenario + replicate + sample, composition_truth, sum)$cell_fraction - 1) > 1e-10) ||
    any(abs(aggregate(rna_contribution_fraction ~ scenario + replicate + sample, composition_truth, sum)$rna_contribution_fraction - 1) > 1e-10)) {
  stop("Cell or RNA contribution fractions do not sum to one")
}
write_gz_csv(data.frame(gene = gene_ids, mixture_y, check.names = FALSE),
             file.path(ANALYSIS, "simulation_mixture_log2_expression.csv.gz"))
write_csv(sample_meta, file.path(ANALYSIS, "simulation_mixture_sample_metadata.csv"))
write_gz_csv(composition_truth, file.path(ANALYSIS, "simulation_mixture_fraction_truth.csv.gz"))

# Algebraic load-bearing truth check: same p, no noise, no post-plant normalization.
p_test <- as.numeric(composition_truth[composition_truth$sample == sample_meta$sample[1], ]$cell_fraction)
linear_control <- drop(profile %*% p_test)
linear_tumor <- drop(profile %*% p_test) * 2^unname(beta)
algebra_check <- data.frame(
  gene = gene_ids,
  beta = unname(beta),
  noiseless_fixed_composition_log2_contrast = log2(linear_tumor / linear_control),
  absolute_error = abs(log2(linear_tumor / linear_control) - unname(beta)),
  stringsAsFactors = FALSE
)
if (max(algebra_check$absolute_error) > 1e-12) stop("Noiseless fixed-composition contrast does not equal planted beta")
write_csv(algebra_check, file.path(ANALYSIS, "simulation_noiseless_truth_identity_check.csv"))

manifest_paths <- c(
  MEAN_PATH, bulk_COMP_PATH, PROTOCOL_MD, PROTOCOL_JSON, SCRIPT_PATH,
  file.path(ANALYSIS, "simulation_gene_truth_and_marker_allocation.csv"),
  file.path(ANALYSIS, "simulation_baseline_linear_profiles.csv.gz"),
  file.path(ANALYSIS, "simulation_synthetic_reference_counts.csv.gz"),
  file.path(ANALYSIS, "simulation_synthetic_reference_cell_metadata.csv"),
  file.path(ANALYSIS, "simulation_mixture_log2_expression.csv.gz"),
  file.path(ANALYSIS, "simulation_mixture_sample_metadata.csv"),
  file.path(ANALYSIS, "simulation_mixture_fraction_truth.csv.gz"),
  file.path(ANALYSIS, "simulation_noiseless_truth_identity_check.csv"),
  file.path(ANALYSIS, "simulation_scenario_centers_and_bulk_contrast.csv")
)
roles <- c("source_background", "scenario_design_reference", "frozen_protocol_md", "frozen_protocol_json", "generator_script",
           rep("frozen_simulation_input", 9L))
manifest <- data.frame(
  path = normalizePath(manifest_paths, winslash = "/", mustWork = TRUE),
  role = roles,
  bytes = file.info(manifest_paths)$size,
  sha256 = vapply(manifest_paths, sha256, character(1)),
  stringsAsFactors = FALSE
)
write_csv(manifest, file.path(ANALYSIS, "simulation_FROZEN_INPUT_MANIFEST.csv"))

cat("simulation inputs frozen before estimator fits.\n")
cat("genes=", nrow(gene_truth), "; markers=", sum(gene_truth$reserved_marker),
    "; target classes=", paste(names(table(gene_truth$gene_class))[names(table(gene_truth$gene_class)) %in% c("null", "composition_only", "expression_only", "both")],
                                as.integer(table(gene_truth$gene_class)[c("null", "composition_only", "expression_only", "both")]), sep = ":", collapse = ","),
    "; mixture samples=", ncol(mixture_y), "; mixture attempts=40.\n", sep = "")
cat("Noiseless identity max absolute error=", format(max(algebra_check$absolute_error), scientific = TRUE), "\n", sep = "")
cat("Separation center max absolute difference=",
    format(max(abs(scenario_defs$separation$control_center - scenario_defs$separation$tumor_center)), digits = 7),
    "; realized mean cell-fraction max absolute difference=",
    format(max(abs(sep_truth_delta)), digits = 7), "\n", sep = "")
cat("Frozen manifest=", file.path(ANALYSIS, "simulation_FROZEN_INPUT_MANIFEST.csv"), "\n", sep = "")
