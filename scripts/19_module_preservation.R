#!/usr/bin/env Rscript
# 19_module_preservation.R
# WGCNA module preservation of discovery (GSE39645) modules in the external
# bulk cohort (GSE108524). Addresses reviewer request for a standard cross-cohort
# module-preservation step.
# Inputs:
#   results/wgcna/wgcna_obj.rds   (list: datExpr [samples x genes], moduleColors)
#   data/processed/GSE108524_expr_gene.csv  (genes x samples)
# Outputs:
#   results/wgcna/module_preservation.csv
#   results/wgcna/fig_module_preservation.png
#   results/wgcna/19_run.log
suppressMessages({library(WGCNA)})
options(stringsAsFactors = FALSE)
set.seed(1)
enableWGCNAThreads()

root <- Sys.getenv("PROJ_ROOT", getwd())
log  <- file(file.path(root, "results/wgcna/19_run.log"), open = "wt")
sink(log, type = "output"); sink(log, type = "message")
cat("== module preservation start:", format(Sys.time()), "==\n")

obj <- readRDS(file.path(root, "results/wgcna/wgcna_obj.rds"))
ref <- obj$datExpr                       # samples x genes (discovery)
modColors <- obj$moduleColors            # length = ncol(ref)
names(modColors) <- colnames(ref)

val <- read.csv(file.path(root, "data/processed/GSE108524_expr_gene.csv"),
                row.names = 1, check.names = FALSE)
val <- as.data.frame(t(val))             # samples x genes

common <- intersect(colnames(ref), colnames(val))
cat("common genes:", length(common), "of", ncol(ref), "discovery module genes\n")
ref <- ref[, common, drop = FALSE]
val <- val[, common, drop = FALSE]
modColors <- modColors[common]

multiExpr  <- list(disc = list(data = ref), valid = list(data = val))
multiColor <- list(disc = modColors)

cat("running modulePreservation (signed, 200 perms)...\n")
mp <- modulePreservation(multiExpr, multiColor,
                         referenceNetworks = 1,
                         networkType = "signed",
                         nPermutations = 200,
                         randomSeed = 1,
                         quickCor = 0,
                         verbose = 3)

stats <- mp$preservation$Z$ref.disc$inColumnsAlsoPresentIn.valid
obsv  <- mp$preservation$observed$ref.disc$inColumnsAlsoPresentIn.valid
out <- data.frame(module    = rownames(stats),
                  moduleSize = obsv$moduleSize,
                  medianRank = mp$preservation$observed$ref.disc$inColumnsAlsoPresentIn.valid$medianRank.pres,
                  Zsummary   = stats$Zsummary.pres,
                  Zdensity   = stats$Zdensity.pres,
                  Zconnectivity = stats$Zconnectivity.pres)
out <- out[order(-out$Zsummary), ]
write.csv(out, file.path(root, "results/wgcna/module_preservation.csv"), row.names = FALSE)
cat("\n== preservation summary ==\n"); print(out)

# Plot Zsummary vs module size (exclude gold/grey reference modules)
keep <- !rownames(stats) %in% c("grey", "gold")
png(file.path(root, "results/wgcna/fig_module_preservation.png"),
    width = 1500, height = 750, res = 150)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
for (i in 1:2) {
  metric <- if (i == 1) out$medianRank else out$Zsummary
  ylab   <- if (i == 1) "Preservation median rank" else "Preservation Zsummary"
  k <- out$module != "grey" & out$module != "gold"
  plot(out$moduleSize[k], metric[k], col = out$module[k], pch = 19, cex = 2,
       xlab = "Module size", ylab = ylab, main = ylab,
       ylim = if (i == 2) range(c(0, metric[k], 12)) else NULL)
  text(out$moduleSize[k], metric[k], labels = out$module[k], pos = 3, cex = 0.8)
  if (i == 2) abline(h = c(2, 10), col = c("blue", "darkgreen"), lty = 2)
}
dev.off()
cat("\n== done:", format(Sys.time()), "==\n")
sink(type = "message"); sink(type = "output")
