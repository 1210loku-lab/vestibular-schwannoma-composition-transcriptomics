#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# N1. Permutation null for cross-cohort (GSE39645 -> GSE108524) reproducibility.
#
# The manuscript reports discovery-vs-validation hub-gene logFC direction
# concordance = 0.9335 and Pearson r = 0.883, with a binomial p that assumes
# independent 50/50 signs. A gene-label permutation null is stronger: it breaks
# the discovery<->validation gene pairing while preserving each cohort's marginal
# logFC distribution (and therefore the sign imbalance), so the null concordance
# reflects agreement expected by chance given the data, not an idealized 0.5.
#
# Independent recomputation: this script recomputes the observed statistics from
# results/validation/hub_logFC_concordance.csv and compares to the stored
# validation_summary.csv values; any mismatch is flagged.
# Outputs only into revision_strengthening_20260818/. Nothing else is touched.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages(library(data.table))
set.seed(20260818)

root   <- Sys.getenv("PROJ_ROOT", getwd())
outdir <- file.path(root, "revision_strengthening_20260818")

conc <- fread(file.path(root, "results/validation/hub_logFC_concordance.csv"))
vsum <- fread(file.path(root, "results/validation/validation_summary.csv"))
stopifnot(all(c("discovery_logFC", "validation_logFC") %in% names(conc)))

d <- conc$discovery_logFC
v <- conc$validation_logFC
n <- length(d)

obs_conc <- mean(sign(d) == sign(v))
obs_r    <- cor(d, v, method = "pearson")

# cross-check against stored summary
stored_conc <- as.numeric(vsum[metric == "logFC_direction_concordance", value])
stored_r    <- as.numeric(vsum[metric == "logFC_pearson_r", value])
concordance_matches <- isTRUE(all.equal(obs_conc, stored_conc, tolerance = 1e-6))
r_matches           <- isTRUE(all.equal(obs_r,    stored_r,    tolerance = 1e-6))
if (!concordance_matches || !r_matches) {
  cat(sprintf("WARN independent recompute mismatch: conc obs=%.6f stored=%.6f | r obs=%.6f stored=%.6f\n",
              obs_conc, stored_conc, obs_r, stored_r))
}

B <- 10000L
null_conc <- numeric(B)
null_r    <- numeric(B)
for (b in seq_len(B)) {
  vp <- v[sample.int(n)]           # permute validation logFC across genes
  null_conc[b] <- mean(sign(d) == sign(vp))
  null_r[b]    <- cor(d, vp)
}

emp_p <- function(null, obs) (1 + sum(null >= obs)) / (length(null) + 1)
p_conc <- emp_p(null_conc, obs_conc)
p_r    <- emp_p(null_r,    obs_r)

summary_tab <- data.table(
  metric = c("n_hub_genes",
             "obs_direction_concordance", "null_concordance_mean", "null_concordance_p95",
             "concordance_perm_p", "concordance_z",
             "obs_pearson_r", "null_r_mean", "null_r_p95", "pearson_r_perm_p", "pearson_r_z",
             "independent_recompute_concordance_matches_stored",
             "independent_recompute_r_matches_stored",
             "n_permutations"),
  value = c(n,
            round(obs_conc, 4), round(mean(null_conc), 4), round(quantile(null_conc, 0.95), 4),
            signif(p_conc, 3), round((obs_conc - mean(null_conc)) / sd(null_conc), 2),
            round(obs_r, 4), round(mean(null_r), 4), round(quantile(null_r, 0.95), 4),
            signif(p_r, 3), round((obs_r - mean(null_r)) / sd(null_r), 2),
            as.integer(concordance_matches), as.integer(r_matches),
            B)
)
fwrite(summary_tab, file.path(outdir, "N1_permutation_null_summary.csv"))

cairo_pdf(file.path(outdir, "N1_permutation_null.pdf"), width = 10.7, height = 5)
par(mfrow = c(1, 2), family = "Arial")
hist(null_conc, breaks = 40, col = "grey80", border = "white",
     xlim = range(c(null_conc, obs_conc)),
     main = "Direction concordance null", xlab = "Concordance (label-permuted)")
abline(v = obs_conc, col = "#c00000", lwd = 3)
legend("topleft", bty = "n",
       legend = sprintf("observed = %.3f\nnull mean = %.3f\nempirical p = 1/%s",
                        obs_conc, mean(null_conc), format(B + 1L, big.mark = ",")))
hist(null_r, breaks = 40, col = "grey80", border = "white",
     xlim = range(c(null_r, obs_r)),
     main = "Pearson r null", xlab = "logFC Pearson r (label-permuted)")
abline(v = obs_r, col = "#c00000", lwd = 3)
legend("topleft", bty = "n",
       legend = sprintf("observed = %.3f\nnull mean = %.3f\nempirical p = 1/%s",
                        obs_r, mean(null_r), format(B + 1L, big.mark = ",")))
dev.off()

cat("=== N1 permutation null ===\n"); print(summary_tab); cat("DONE N1\n")
