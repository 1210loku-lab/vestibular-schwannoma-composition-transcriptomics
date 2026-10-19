#!/usr/bin/env Rscript

ROOT <- Sys.getenv("VS_INPUT_ROOT", unset = ".")
EXPORT_ROOT <- Sys.getenv("VS_OUTPUT_ROOT", unset = ".")
ANALYSIS <- file.path(EXPORT_ROOT, "outputs")
FIGURES <- file.path(EXPORT_ROOT, "figures")
suppressPackageStartupMessages({
  library(ggplot2)
})

source_data <- read.csv(file.path(ANALYSIS, "donor_Figure_S17_source.csv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
panel_order <- c("NRG1", "L1CAM", "FCGBP", "PRRX1")
group_levels <- c("Nerve control", "Vestibular schwannoma (VS)")
group_palette <- c("Nerve control" = "#a6bf9c", "Vestibular schwannoma (VS)" = "#9fb3d8")
source_data$group <- factor(source_data$group, levels = group_levels)

panel_plots <- lapply(panel_order, function(gene) {
  dat <- source_data[source_data$panel_gene == gene, , drop = FALSE]
  if (!nrow(dat) || length(unique(dat$donor_id)) != 9L) {
    stop("Figure source is incomplete for predeclared panel: ", gene)
  }
  lineage <- dat$panel_lineage[1]
  summary <- do.call(rbind, lapply(group_levels, function(g) {
    vals <- dat$log2_normalized_count_plus_1[dat$group == g]
    data.frame(group = factor(g, levels = group_levels), mean = mean(vals), SD = sd(vals))
  }))
  ggplot(dat, aes(x = group, y = log2_normalized_count_plus_1, color = group)) +
    geom_col(data = summary[summary$group == "Nerve control", , drop = FALSE],
      aes(x = group, y = mean), inherit.aes = FALSE, width = 0.52,
      fill = scales::alpha(group_palette[["Nerve control"]], 0.32),
      color = group_palette[["Nerve control"]], linewidth = 0.45) +
    geom_col(data = summary[summary$group == "Vestibular schwannoma (VS)", , drop = FALSE],
      aes(x = group, y = mean), inherit.aes = FALSE, width = 0.52,
      fill = scales::alpha(group_palette[["Vestibular schwannoma (VS)"]], 0.32),
      color = group_palette[["Vestibular schwannoma (VS)"]], linewidth = 0.45) +
    geom_errorbar(data = summary,
      aes(x = group, ymin = mean - SD, ymax = mean + SD, color = group),
      inherit.aes = FALSE, width = 0.14, linewidth = 0.45) +
    geom_point(position = position_jitter(width = 0.075, height = 0, seed = 20261015),
      size = 1.9, alpha = 0.95, stroke = 0) +
    scale_color_manual(values = group_palette, drop = FALSE) +
    scale_x_discrete(labels = c("Nerve control" = "Nerve control\nn=2",
      "Vestibular schwannoma (VS)" = "VS\nn=7")) +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.06))) +
    labs(
      title = bquote(italic(.(gene)) ~ "—" ~ .(lineage)),
      x = NULL,
      y = "log2(DESeq2-normalized count + 1)"
    ) +
    theme_classic(base_size = 8, base_family = "Arial") +
    theme(
      axis.text = element_text(family = "Arial", size = 7, color = "#222222"),
      axis.title.y = element_text(family = "Arial", size = 7),
      plot.title = element_text(family = "Arial", size = 8, face = "plain", hjust = 0),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(4, 5, 4, 5)
    )
})

s17 <- patchwork::wrap_plots(panel_plots, ncol = 2) +
  patchwork::plot_annotation(
    tag_levels = "a",
    title = "Donor-level candidate expression in broad compartments",
    subtitle = "Points are biological donors; bars are group means and whiskers are SD. Values are DESeq2-normalized counts on a log2(count + 1) scale.",
    caption = paste0(
      "Panels fixed a priori: NRG1/L1CAM (Schwann); FCGBP (myeloid);\n",
      "PRRX1 (fibroblast). Schwann combines normal and tumor-associated states;\n",
      "this is not a pure neoplastic-cell compartment. Panel-specific y-axis ranges."
    ),
    theme = theme(
      plot.title = element_text(family = "Arial", face = "bold", size = 10),
      plot.subtitle = element_text(family = "Arial", size = 7),
      plot.caption = element_text(family = "Arial", size = 6, hjust = 0,
                                  lineheight = 1.05, margin = margin(t = 4)),
      plot.tag = element_text(family = "Arial", face = "bold", size = 10)
    )
  )

stem <- file.path(FIGURES, "Figure_S17_legacy_donor_pseudobulk")
svglite::svglite(paste0(stem, ".svg"), width = 7.5, height = 6.2)
print(s17)
grDevices::dev.off()
grDevices::cairo_pdf(paste0(stem, ".pdf"), width = 7.5, height = 6.2, family = "Arial")
print(s17)
grDevices::dev.off()
ragg::agg_png(paste0(stem, ".png"), width = 7.5, height = 6.2, units = "in", res = 600)
print(s17)
grDevices::dev.off()
message("S17 refreshed from saved donor-expression table: ", stem)
