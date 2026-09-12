#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

project_root <- normalizePath(Sys.getenv("PROJ_ROOT", getwd()), mustWork = TRUE)
source_dir <- file.path(project_root, "results", "deconv_multimethod")
crossref_file <- file.path(
  project_root,
  "results", "robustness", "cross_reference",
  "global_compartment_concordance_summary.csv"
)
output_dir <- Sys.getenv(
  "CBC_FIGURE_DIR",
  file.path(project_root, "submission", "CBC_20260911_final", "figures", "main")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stability <- read.csv(
  file.path(source_dir, "reference_by_method_stability.csv"),
  check.names = FALSE
)
calibration <- read.csv(
  file.path(source_dir, "known_mixture_calibration.csv"),
  check.names = FALSE
)
rare <- read.csv(
  file.path(source_dir, "known_mixture_rare_lowprop.csv"),
  check.names = FALSE
)
crossref <- read.csv(crossref_file, check.names = FALSE)

cell_order <- c(
  "Schwann/tumour", "Myeloid", "T cell",
  "Fibroblast", "Endothelial", "Pericyte"
)
cell_labels <- c(
  "Schwann/tumour" = "Schwann/tumour",
  "Myeloid" = "Myeloid",
  "T cell" = "T cell",
  "Fibroblast" = "Fibroblast",
  "Endothelial" = "Endothelial",
  "Pericyte" = "Pericyte"
)

theme_cbc <- theme_classic(base_family = "Arial", base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0),
    plot.subtitle = element_text(size = 8.5, colour = "#333333"),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 8, colour = "black"),
    strip.background = element_rect(fill = "#F2F2F2", colour = "#BDBDBD"),
    strip.text = element_text(face = "bold", size = 8),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 8),
    plot.margin = margin(7, 8, 7, 8)
  )

# Panel A: keep method and reference effects explicit and give every label room.
a_method <- stability[c("cohort", "cell_type", "method_effect_med")]
names(a_method)[3] <- "rho"
a_method$effect <- "Method effect"
a_ref <- stability[c("cohort", "cell_type", "reference_effect_med")]
names(a_ref)[3] <- "rho"
a_ref$effect <- "Reference effect"
a_data <- rbind(a_method, a_ref)
a_data$cell_type <- factor(a_data$cell_type, levels = rev(cell_order))
a_data$configuration <- factor(
  paste(a_data$cohort, a_data$effect, sep = "\n"),
  levels = c(
    "GSE39645\nMethod effect", "GSE39645\nReference effect",
    "GSE108524\nMethod effect", "GSE108524\nReference effect"
  )
)

p_a <- ggplot(a_data, aes(configuration, cell_type, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f", rho)), family = "Arial", size = 3) +
  scale_fill_gradientn(
    colours = c("#B2182B", "#FDD49E", "#D9F0D3", "#1B7837"),
    limits = c(-0.1, 1), name = "Spearman rho"
  ) +
  scale_y_discrete(labels = cell_labels) +
  labs(
    tag = "A",
    title = "Cross-configuration rank concordance",
    subtitle = "Median per-sample Spearman correlations within cohort and compartment",
    x = NULL, y = NULL
  ) +
  theme_cbc +
  theme(axis.text.x = element_text(size = 8, lineheight = 0.9))

# Panel B: display every metric named in the legend.
b_ccc <- calibration[c("method", "cell_type", "CCC")]
names(b_ccc)[3] <- "value"
b_ccc$metric <- "CCC"
b_mae <- calibration[c("method", "cell_type", "MAE")]
names(b_mae)[3] <- "value"
b_mae$value <- 100 * b_mae$value
b_mae$metric <- "Overall MAE (%)"
b_rare <- rare[c("method", "cell_type", "MAE_rare")]
names(b_rare)[3] <- "value"
b_rare$value <- 100 * b_rare$value
b_rare$metric <- "Rare-proportion MAE (%)"
b_data <- rbind(b_ccc, b_mae, b_rare)
b_data$cell_type <- factor(b_data$cell_type, levels = cell_order)
b_data$method <- factor(b_data$method, levels = c("NNLS", "dtangle"))
b_data$metric <- factor(
  b_data$metric,
  levels = c("CCC", "Overall MAE (%)", "Rare-proportion MAE (%)")
)

p_b <- ggplot(b_data, aes(cell_type, value, fill = method)) +
  geom_col(position = position_dodge(width = 0.76), width = 0.68) +
  facet_wrap(~metric, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c("NNLS" = "#3B6FB6", "dtangle" = "#D6604D")) +
  scale_x_discrete(labels = cell_labels) +
  labs(
    tag = "B",
    title = "Donor-held-out synthetic-mixture recovery",
    subtitle = "GSE216783; 13 donors and 100 mixtures per held-out donor",
    x = NULL, y = NULL
  ) +
  theme_cbc +
  theme(
    axis.text.x = element_text(angle = 38, hjust = 1, vjust = 1, size = 7),
    panel.spacing.x = unit(10, "pt")
  )

# Panel C: use the values directly from the cross-reference source table.
c_data <- crossref[crossref$min_top_mean == 0, ]
c_data <- c_data[match(c(1, 1.2, 1.5, 2), c_data$min_dominance_ratio), ]
c_data$threshold <- factor(
  c("All", "≥1.2", "≥1.5", "≥2.0"),
  levels = c("All", "≥1.2", "≥1.5", "≥2.0")
)
c_data$label <- sprintf("%.1f%%\n(n=%d)", 100 * c_data$concordance, c_data$n_evaluable)

p_c <- ggplot(c_data, aes(threshold, 100 * concordance, group = 1)) +
  geom_line(colour = "#1B7837", linewidth = 0.8) +
  geom_point(shape = 21, size = 3.2, stroke = 0.8, fill = "white", colour = "#1B7837") +
  geom_text(aes(label = label), family = "Arial", size = 2.8, vjust = -0.8) +
  scale_y_continuous(limits = c(55, 108), breaks = c(60, 70, 80, 90, 100)) +
  labs(
    tag = "C",
    title = "Cross-reference assignment agreement",
    subtitle = "Agreement rises with top-to-second dominance",
    x = "Top-to-second expression ratio", y = "Agreement (%)"
  ) +
  theme_cbc

combined <- (p_a | p_c) / p_b +
  plot_layout(heights = c(1, 1.15), widths = c(1.45, 1)) &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 13))

pdf_file <- file.path(output_dir, "Fig5.pdf")
svg_file <- file.path(output_dir, "Fig5_editable.svg")
ggsave(pdf_file, combined, width = 12, height = 7.2, units = "in", device = cairo_pdf)
ggsave(svg_file, combined, width = 12, height = 7.2, units = "in", device = svglite::svglite)

message("Created: ", pdf_file)
message("Created: ", svg_file)
