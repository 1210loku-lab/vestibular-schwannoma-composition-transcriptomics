#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# N3. Evidence-tiered candidate resource.
#
# Turns the 354-transcript cross-evidence map into an explicit, machine-readable
# tiered resource with transparent, reproducible rules. Tiers rank *evidence
# strength for follow-up*, NOT effect size, mechanism, or causal role.
#
# Five evidence layers (already in the map): ev_bulk_sig, ev_hub, ev_cross_conc
# (cross-cohort direction concordance), ev_pb_intrinsic (within-compartment
# pseudobulk signal), ev_adj_retain (survives composition adjustment).
# evidence_layer_count = their sum (2-5 across the candidate set).
#
# Tier rules (composition-aware layer is the key discriminator):
#   Tier 1  : all 5 layers (bulk-sig + hub + cross-cohort concordant +
#             within-compartment intrinsic + composition-retained)
#   Tier 2  : 4 layers AND cross-cohort concordant
#   Tier 3  : remaining concordant/hub candidates (2-3 layers; typically
#             composition-attenuated or lacking intrinsic-compartment support)
#
# Robustness supports (leave-one-VS, composition-multiverse joint threshold,
# dual single-cell cross-reference) are attached as within-tier annotations.
# Independent check: recomputes evidence_layer_count from the five flags and
# compares to the stored column. Outputs only into revision_strengthening_20260818/.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
options(stringsAsFactors = FALSE)

root   <- Sys.getenv("PROJ_ROOT", getwd())
outdir <- file.path(root, "revision_strengthening_20260818")
map    <- fread(file.path(root, "results/robustness/TableS1_compartment_cross_evidence_map_robustness.csv"))

layers <- c("ev_bulk_sig", "ev_hub", "ev_cross_conc", "ev_pb_intrinsic", "ev_adj_retain")
map[, recomputed_layer_count := rowSums(as.data.table(lapply(.SD, function(x) as.integer(x == TRUE))),
                                        na.rm = TRUE), .SDcols = layers]
layercount_reproduced <- all(map$recomputed_layer_count == map$evidence_layer_count, na.rm = TRUE)
if (!layercount_reproduced) cat("WARN: recomputed evidence_layer_count != stored column\n")

# Tier labels are neutral (layer-count based); layers are partly analytically
# dependent (hub membership includes the DEG criterion), so tiers are descriptive
# prioritization criteria, not independent-evidence accumulation.
map[, evidence_tier := fifelse(
      recomputed_layer_count >= 5, "Tier1_five_layer",
      fifelse(recomputed_layer_count == 4 & cross_cohort_concordant == TRUE,
              "Tier2_four_layer_concordant", "Tier3_exploratory"))]

# Hard validations (fail loudly rather than silently pass)
stopifnot(
  nrow(map) == 354L,
  anyDuplicated(map$gene) == 0L,
  !anyNA(map$evidence_tier),
  identical(as.logical(map$ev_cross_conc), as.logical(map$cross_cohort_concordant))
)

# Robustness annotation flags (within-tier, do not change tier)
map[, robust_leave_one_vs := leave_one_vs_direction_support == 1]
map[, robust_composition_multiverse := baseline_joint_threshold_support >= 0.8]
map[, robust_dual_singlecell := cross_reference_concordant == TRUE]

resource <- map[, .(
  gene,
  evidence_tier,
  evidence_layer_count = recomputed_layer_count,
  ev_bulk_sig, ev_hub, ev_cross_conc, ev_pb_intrinsic, ev_adj_retain,
  bulk_logFC, cross_cohort_concordant,
  dominant_compartment = pb_dominant_celltype,
  composition_adj_class,
  robust_leave_one_vs, robust_composition_multiverse, robust_dual_singlecell
)][order(factor(evidence_tier, levels = c("Tier1_five_layer","Tier2_four_layer_concordant","Tier3_exploratory")),
         -evidence_layer_count, -abs(bulk_logFC))]

fwrite(resource, file.path(outdir, "N3_candidate_evidence_tier_resource.csv"))

tier_summary <- resource[, .(
  n = .N,
  n_composition_retained = sum(ev_adj_retain == TRUE, na.rm = TRUE),
  n_within_compartment_intrinsic = sum(ev_pb_intrinsic == TRUE, na.rm = TRUE),
  n_cross_cohort_concordant = sum(cross_cohort_concordant == TRUE, na.rm = TRUE),
  n_robust_all3 = sum(robust_leave_one_vs & robust_composition_multiverse & robust_dual_singlecell, na.rm = TRUE)
), by = evidence_tier][order(evidence_tier)]
fwrite(tier_summary, file.path(outdir, "N3_tier_summary.csv"))

# Named Tier 1 candidates for the manuscript resource callout
tier1 <- resource[evidence_tier == "Tier1_five_layer",
                  .(gene, dominant_compartment, bulk_logFC = round(bulk_logFC, 2),
                    robust_all3 = robust_leave_one_vs & robust_composition_multiverse & robust_dual_singlecell)]
fwrite(tier1, file.path(outdir, "N3_tier1_candidates.csv"))

cat("=== N3 evidence-tier resource ===\n")
cat("layer_count independently reproduced stored column:", layercount_reproduced, "\n\n")
print(tier_summary)
cat("\n-- Tier 1 candidates --\n"); print(tier1)
cat("DONE N3\n")
