#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(WGCNA)
})

options(stringsAsFactors = FALSE)
set.seed(42)
allowWGCNAThreads()
cor <- WGCNA::cor

root <- Sys.getenv("PROJ_ROOT", getwd())
outdir <- file.path(root, "results/hub")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

obj <- readRDS(file.path(root, "results/wgcna/wgcna_obj.rds"))
deg <- fread(file.path(root, "results/deg/deg_full_gene.csv"))

stopifnot(
  all(c("datExpr", "moduleColors", "trait") %in% names(obj)),
  ncol(obj$datExpr) == length(obj$moduleColors),
  all(c("gene", "logFC", "adj.P.Val") %in% names(deg))
)

datExpr <- as.data.frame(obj$datExpr)
moduleColors <- as.character(obj$moduleColors)
names(moduleColors) <- colnames(datExpr)
trait <- as.numeric(obj$trait)

MEs <- orderMEs(moduleEigengenes(datExpr, colors = moduleColors)$eigengenes)
kme <- signedKME(datExpr, MEs, outputColumnName = "kME")
colnames(kme) <- sub("^kME", "", colnames(kme))
gs <- as.numeric(cor(datExpr, trait, use = "pairwise.complete.obs"))
names(gs) <- colnames(datExpr)

module_specs <- data.table(
  module = c("turquoise", "yellow"),
  direction = c("up", "down")
)

all_module_genes <- rbindlist(lapply(seq_len(nrow(module_specs)), function(i) {
  mod <- module_specs$module[i]
  direction <- module_specs$direction[i]
  genes <- colnames(datExpr)[moduleColors == mod]
  stopifnot(mod %in% colnames(kme), length(genes) > 0)

  tab <- data.table(
    gene = genes,
    module = mod,
    direction = direction,
    kME = kme[genes, mod],
    GS = gs[genes]
  )
  tab <- merge(
    tab,
    deg[, .(gene, logFC, adj.P.Val)],
    by = "gene",
    all.x = TRUE,
    sort = FALSE
  )
  tab[, is_hub := (
    abs(kME) > 0.8 &
      abs(GS) > 0.5 &
      !is.na(adj.P.Val) & adj.P.Val < 0.05 &
      !is.na(logFC) & abs(logFC) > 1
  )]
  tab <- tab[order(-abs(kME), -abs(GS))]

  p <- ggplot(tab, aes(x = kME, y = GS, color = is_hub)) +
    geom_hline(yintercept = c(-0.5, 0.5), linetype = 2, color = "grey60") +
    geom_vline(xintercept = c(-0.8, 0.8), linetype = 2, color = "grey60") +
    geom_point(alpha = 0.75, size = 1.8) +
    scale_color_manual(values = c(`FALSE` = "grey75", `TRUE` = "#D55E00")) +
    labs(
      title = sprintf("%s module: kME versus gene significance", mod),
      x = sprintf("signed kME (%s)", mod),
      y = "Gene significance (correlation with Tumor)",
      color = "Hub criterion"
    ) +
    theme_bw(base_size = 11)
  ggsave(
    file.path(outdir, sprintf("fig_kME_GS_%s.png", mod)),
    p,
    width = 6.5,
    height = 5.2,
    dpi = 180
  )
  tab
}), use.names = TRUE)

hubs <- all_module_genes[is_hub == TRUE, .(
  gene, module, direction, kME, GS, logFC, adj.P.Val
)]
hubs <- hubs[order(module, -abs(kME), -abs(GS))]

hub_up <- hubs[module == "turquoise"]
hub_down <- hubs[module == "yellow"]

fwrite(hub_up, file.path(outdir, "hub_up_turquoise.csv"))
fwrite(hub_down, file.path(outdir, "hub_down_yellow.csv"))
fwrite(hubs, file.path(outdir, "hub_all.csv"))
fwrite(
  all_module_genes[, .(
    gene, module, direction, kME, GS, logFC, adj.P.Val, is_hub
  )],
  file.path(outdir, "module_gene_metrics.csv")
)

summary_tab <- rbindlist(lapply(c("turquoise", "yellow"), function(mod) {
  x <- hubs[module == mod]
  data.table(
    module = mod,
    direction = unique(module_specs[module == mod, direction]),
    hub_count = nrow(x),
    top10_genes = paste(head(x$gene, 10), collapse = ";")
  )
}))
fwrite(summary_tab, file.path(outdir, "hub_summary.csv"))

print(summary_tab)
cat("DONE 05\n")
