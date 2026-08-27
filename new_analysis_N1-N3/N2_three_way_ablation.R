#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# N2. Three-way ablation: how candidate ranking and biological interpretation
# change under three increasingly composition-aware strategies.
#
#   (A) Naive bulk DE            : original_logFC / original_adj.P.Val
#   (B) Composition-adjusted DE  : adjusted_logFC / adjusted_adj.P.Val
#                                  (cell-fraction covariates from deconvolution)
#   (C) Compartment-aware        : single-cell dominant-compartment assignment
#                                  (pb_dominant_celltype in the cross-evidence map)
#
# Message being quantified (not a novelty claim, not causal): naive bulk ranking
# and composition-adjusted ranking diverge for a large block of transcripts.
# Here, "attenuated" is a machine-readable legacy label meaning only that a
# transcript lost the prespecified significance thresholds after adjustment;
# it does not identify a composition-driven effect. Compartment assignment is
# an annotation layer rather than a third inferential model.
#
# Independent check: the A->B classification counts are recomputed here from the
# raw logFC/adjP and compared to results/deconv/deg_adjustment_summary.csv.
# Outputs only into revision_strengthening_20260818/.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
options(stringsAsFactors = FALSE)

root   <- Sys.getenv("PROJ_ROOT", getwd())
outdir <- file.path(root, "revision_strengthening_20260818")

adj  <- fread(file.path(root, "results/deconv/deg_composition_adjusted.csv"))
smry <- fread(file.path(root, "results/deconv/deg_adjustment_summary.csv"))
map  <- fread(file.path(root, "results/robustness/TableS1_compartment_cross_evidence_map_robustness.csv"))

SIG_P  <- 0.05
LFC_AB <- 1

# --- Independent recompute of A vs B significance classification ------------
adj[, A_sig := original_adj.P.Val < SIG_P & abs(original_logFC) > LFC_AB]
adj[, B_sig := adjusted_adj.P.Val < SIG_P & abs(adjusted_logFC) > LFC_AB]
adj[, recomputed_class := fifelse( A_sig &  B_sig, "retained",
                          fifelse( A_sig & !B_sig, "attenuated",
                          fifelse(!A_sig &  B_sig, "adjusted-only", "neither")))]
recomp <- adj[, .N, by = recomputed_class][order(recomputed_class)]
setnames(recomp, "N", "recomputed_Freq")
cmp <- merge(smry, recomp, by.x = "classification", by.y = "recomputed_class", all = TRUE)
cmp[, matches_stored := Freq == recomputed_Freq]
class_reproduced <- all(cmp$matches_stored, na.rm = TRUE)
if (!class_reproduced) cat("WARN: A->B classification did not reproduce stored deg_adjustment_summary.csv\n")

# --- Ranking divergence A vs B ---------------------------------------------
adj[, rankmetric_A := -log10(pmax(original_adj.P.Val, 1e-300))]
adj[, rankmetric_B := -log10(pmax(adjusted_adj.P.Val, 1e-300))]
sp_all  <- cor(adj$rankmetric_A, adj$rankmetric_B, method = "spearman")
naiveA  <- adj[A_sig == TRUE]
sp_sig  <- cor(naiveA$rankmetric_A, naiveA$rankmetric_B, method = "spearman")
lfc_sp  <- cor(adj$original_logFC, adj$adjusted_logFC, method = "spearman")

# --- Layer C: compartment reassignment of the A->B classes -----------------
mm <- merge(adj[, .(gene, A_sig, B_sig, recomputed_class)],
            map[, .(gene, pb_dominant_celltype)], by = "gene")
comp_tab <- mm[recomputed_class %in% c("retained", "attenuated"),
               .N, by = .(recomputed_class, pb_dominant_celltype)][order(recomputed_class, -N)]

summary_tab <- data.table(
  metric = c("n_genes_total",
             "A_significant_naive", "B_significant_adjusted",
             "class_retained", "class_attenuated", "class_adjusted_only",
             "spearman_rank_A_vs_B_all", "spearman_rank_A_vs_B_naiveSig",
             "spearman_logFC_A_vs_B",
             "candidates_with_compartment_assignment",
             "AtoB_classification_reproduced_stored"),
  value = c(nrow(adj),
            sum(adj$A_sig), sum(adj$B_sig),
            recomp[recomputed_class == "retained", recomputed_Freq],
            recomp[recomputed_class == "attenuated", recomputed_Freq],
            recomp[recomputed_class == "adjusted-only", recomputed_Freq],
            round(sp_all, 3), round(sp_sig, 3), round(lfc_sp, 3),
            nrow(mm), as.integer(class_reproduced))
)

fwrite(summary_tab, file.path(outdir, "N2_ablation_summary.csv"))
fwrite(cmp,         file.path(outdir, "N2_classification_reproduce_check.csv"))
fwrite(comp_tab,    file.path(outdir, "N2_compartment_by_class.csv"))

# --- Figure: A vs B logFC scatter (coloured by class) + class bar ----------
pal <- c(retained = "#0072B2", attenuated = "#D55E00",
         `adjusted-only` = "#009E73", neither = "grey80")
cairo_pdf(file.path(outdir, "N2_ablation.pdf"), width = 10.7, height = 5.3)
par(mfrow = c(1, 2), family = "Arial", mar = c(5.2, 4.6, 3.6, 1.2))
plot(adj$original_logFC, adj$adjusted_logFC,
     col = pal[adj$recomputed_class], pch = 16, cex = 0.5,
     xlab = "Naive bulk logFC (A)", ylab = "Composition-adjusted logFC (B)",
     main = sprintf("A vs B logFC (Spearman %.2f)", lfc_sp))
abline(0, 1, col = "grey40", lty = 2); abline(h = 0, v = 0, col = "grey85")
legend("topleft", bty = "n", pch = 16, col = pal[c("retained","attenuated","adjusted-only")],
       legend = c("retained (A & B)", "lost support (A only)", "adjusted only (B only)"))
bars <- recomp[recomputed_class != "neither"]
display_names <- c(`adjusted-only` = "adjusted only",
                   attenuated = "lost support", retained = "retained")
barplot(bars$recomputed_Freq, names.arg = display_names[bars$recomputed_class],
        col = pal[bars$recomputed_class], las = 1, cex.names = 0.82,
        main = "A->B significance reclassification", ylab = "n transcripts")
dev.off()

cat("=== N2 three-way ablation ===\n"); print(summary_tab)
cat("\n-- compartment x class (top rows) --\n"); print(head(comp_tab, 12))
cat("DONE N2\n")
