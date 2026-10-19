#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
script_dir <- dirname(script_path)
package_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
input_root <- normalizePath(Sys.getenv("VS_INPUT_ROOT", unset = file.path(package_root, "inputs")), mustWork = TRUE)
output_root <- normalizePath(Sys.getenv("VS_OUTPUT_ROOT", unset = package_root), mustWork = TRUE)
output_dir <- file.path(output_root, "figures")
plotdata_dir <- file.path(output_root, "outputs", "figure_plotted_data")
qa_dir <- file.path(output_root, "outputs", "figure_qa", "figure_3")
source_manifest_dir <- file.path(output_root, "outputs", "figure_source_manifests")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plotdata_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_manifest_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(systemfonts)
})
poppler_tools <- c("pdfinfo", "pdftotext", "pdftoppm")
missing_poppler <- poppler_tools[!nzchar(Sys.which(poppler_tools))]
if (length(missing_poppler)) stop("Missing Poppler utilities on PATH: ", paste(missing_poppler, collapse = ", "))

source_files <- c(
  reference_by_method_stability = file.path("results", "deconv_multimethod", "reference_by_method_stability.csv"),
  known_mixture_calibration = file.path("results", "deconv_multimethod", "known_mixture_calibration.csv"),
  known_mixture_rare_lowprop = file.path("results", "deconv_multimethod", "known_mixture_rare_lowprop.csv"),
  cross_reference_concordance = file.path("results", "robustness", "cross_reference", "global_compartment_concordance_summary.csv")
)
source_paths <- setNames(file.path(input_root, unname(source_files)), names(source_files))
missing_sources <- source_paths[!file.exists(source_paths)]
if (length(missing_sources)) stop("Missing source table(s): ", paste(missing_sources, collapse = "; "))

stability <- read.csv(source_paths[["reference_by_method_stability"]], check.names = FALSE)
calibration <- read.csv(source_paths[["known_mixture_calibration"]], check.names = FALSE)
rare <- read.csv(source_paths[["known_mixture_rare_lowprop"]], check.names = FALSE)
crossref <- read.csv(source_paths[["cross_reference_concordance"]], check.names = FALSE)

cell_order <- c(
  "Schwann/tumour", "Myeloid", "T cell",
  "Fibroblast", "Endothelial", "Pericyte"
)
display_labels <- c(
  "Schwann/tumour" = "Schwann/\ntumor",
  "Myeloid" = "Myeloid",
  "T cell" = "T cell",
  "Fibroblast" = "Fibroblast",
  "Endothelial" = "Endothelial",
  "Pericyte" = "Pericyte"
)
theme_fig3 <- theme_classic(base_family = "Arial", base_size = 8) +
  theme(
    plot.title = element_text(
      family = "Arial", face = "bold", size = 9, colour = "#202020",
      hjust = 0, margin = margin(l = 13, b = 1, unit = "pt")
    ),
    plot.subtitle = element_text(
      family = "Arial", size = 8, colour = "#333333",
      hjust = 0, margin = margin(l = 13, b = 2, unit = "pt")
    ),
    plot.title.position = "plot",
    plot.tag = element_text(family = "Arial", face = "bold", size = 12, colour = "#111111"),
    plot.tag.position = c(0, 1),
    plot.margin = margin(t = 1, r = 3, b = 0.6, l = 3, unit = "mm"),
    axis.title = element_text(family = "Arial", size = 8, colour = "#222222"),
    axis.text = element_text(family = "Arial", size = 8, colour = "#222222"),
    axis.line = element_line(linewidth = 0.35, colour = "#222222"),
    axis.ticks = element_line(linewidth = 0.35, colour = "#222222"),
    legend.title = element_text(family = "Arial", size = 8, colour = "#222222"),
    legend.text = element_text(family = "Arial", size = 8, colour = "#222222"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(7, "mm"),
    strip.background = element_rect(fill = "#F2F2F2", colour = "#BDBDBD", linewidth = 0.35),
    strip.text = element_text(family = "Arial", face = "bold", size = 8, colour = "#222222"),
    panel.grid = element_blank()
  )

# Panel a: medians of the source table's per-sample Spearman rank correlations.
a_method <- stability[c("cohort", "cell_type", "method_effect_med")]
names(a_method)[3] <- "rho"
a_method$effect <- "Method effect"
a_reference <- stability[c("cohort", "cell_type", "reference_effect_med")]
names(a_reference)[3] <- "rho"
a_reference$effect <- "Reference effect"
a_data <- rbind(a_method, a_reference)
a_data$configuration <- factor(
  paste(a_data$cohort, a_data$effect, sep = "\n"),
  levels = c(
    "GSE39645\nMethod effect", "GSE39645\nReference effect",
    "GSE108524\nMethod effect", "GSE108524\nReference effect"
  )
)
a_data$cell_type <- factor(a_data$cell_type, levels = rev(cell_order))
a_data$cell_type_display <- unname(display_labels[as.character(a_data$cell_type)])
a_data$label_colour <- ifelse(a_data$rho < 0.15 | a_data$rho > 0.78, "white", "#111111")
a_data$source_table <- basename(source_files[["reference_by_method_stability"]])

p_a <- ggplot(a_data, aes(configuration, cell_type, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_tile(
    data = subset(a_data, cohort == "GSE108524" & cell_type == "T cell" & effect == "Method effect"),
    inherit.aes = FALSE, aes(configuration, cell_type),
    fill = NA, colour = "#202020", linewidth = 0.55
  ) +
  geom_text(
    aes(label = sprintf("%.2f", rho), colour = label_colour),
    family = "Arial", size = 8 / .pt, show.legend = FALSE
  ) +
  scale_colour_identity() +
  scale_fill_gradientn(
    colours = c("#B2182B", "#FDD49E", "#D9F0D3", "#1B7837"),
    limits = c(-0.1, 1), breaks = c(0, 0.25, 0.5, 0.75, 1),
    oob = scales::squish,
    name = expression(Spearman~rho),
    guide = guide_colourbar(
      direction = "horizontal", title.position = "top", title.hjust = 0.5,
      barwidth = unit(43, "mm"), barheight = unit(2.7, "mm"),
      ticks.colour = "#333333"
    )
  ) +
  scale_x_discrete(expand = expansion(add = 0.02)) +
  scale_y_discrete(labels = display_labels, expand = expansion(add = 0.02)) +
  labs(
    tag = "a",
    title = "Cross-configuration rank concordance",
    subtitle = "Median per-sample Spearman ρ within cohort and compartment",
    x = NULL, y = NULL
  ) +
  theme_fig3 +
  theme(
    axis.text.x = element_text(size = 8, lineheight = 0.92, colour = "#222222"),
    axis.text.y = element_text(size = 8, lineheight = 0.92, colour = "#222222"),
    legend.margin = margin(t = 0, b = 0, unit = "pt"),
    legend.box.margin = margin(t = 0, b = 0, unit = "pt")
  )

# Panel b: direct display of frozen CCC and MAE summaries; MAE is expressed in percent.
b_ccc <- calibration[c("method", "cell_type", "CCC")]
names(b_ccc)[3] <- "value"
b_ccc$metric <- "CCC"
b_ccc$source_table <- basename(source_files[["known_mixture_calibration"]])
b_ccc$n <- NA_integer_
b_mae <- calibration[c("method", "cell_type", "MAE")]
names(b_mae)[3] <- "value"
b_mae$value <- 100 * b_mae$value
b_mae$metric <- "Overall MAE (%)"
b_mae$source_table <- basename(source_files[["known_mixture_calibration"]])
b_mae$n <- NA_integer_
b_rare <- rare[c("method", "cell_type", "MAE_rare", "n")]
names(b_rare)[3] <- "value"
b_rare$value <- 100 * b_rare$value
b_rare$metric <- "Rare-proportion MAE (%)"
b_rare$source_table <- basename(source_files[["known_mixture_rare_lowprop"]])
b_data <- rbind(
  b_ccc[c("method", "cell_type", "value", "metric", "source_table", "n")],
  b_mae[c("method", "cell_type", "value", "metric", "source_table", "n")],
  b_rare[c("method", "cell_type", "value", "metric", "source_table", "n")]
)
b_data$method <- factor(b_data$method, levels = c("NNLS", "dtangle"))
b_data$cell_type <- factor(b_data$cell_type, levels = rev(cell_order))
b_data$cell_type_display <- unname(display_labels[as.character(b_data$cell_type)])
b_data$metric <- factor(
  b_data$metric,
  levels = c("CCC", "Overall MAE (%)", "Rare-proportion MAE (%)")
)
b_data$y_center <- as.numeric(b_data$cell_type)
b_data$y <- b_data$y_center + ifelse(b_data$method == "NNLS", 0.15, -0.15)
b_wide <- reshape(
  b_data[c("cell_type", "metric", "method", "value")],
  idvar = c("cell_type", "metric"), timevar = "method", direction = "wide"
)
b_wide$y_center <- as.numeric(factor(b_wide$cell_type, levels = rev(cell_order)))
b_wide$source_table <- "known_mixture_calibration.csv + known_mixture_rare_lowprop.csv"

p_b <- ggplot(b_data, aes(value, y)) +
  geom_segment(
    data = b_wide,
    inherit.aes = FALSE,
    aes(x = value.NNLS, xend = value.dtangle, y = y_center + 0.15, yend = y_center - 0.15),
    colour = "#A5A5A5", linewidth = 0.4
  ) +
  geom_point(aes(colour = method, shape = method), size = 2.3, stroke = 0.45) +
  facet_wrap(~metric, nrow = 1, scales = "free_x") +
  scale_colour_manual(
    values = c("NNLS" = "#3B6FB6", "dtangle" = "#D6604D"),
    breaks = c("NNLS", "dtangle"), name = NULL
  ) +
  scale_shape_manual(values = c("NNLS" = 16, "dtangle" = 15), name = NULL, guide = "none") +
  scale_y_continuous(
    breaks = seq_along(rev(cell_order)), labels = unname(display_labels[rev(cell_order)]),
    limits = c(0.5, 6.5), expand = expansion(mult = 0)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.08))) +
  labs(
    tag = "b",
    title = "Donor-held-out synthetic-mixture recovery",
    subtitle = "GSE216783; within-dataset holdout: 13 donors × 100 mixtures/donor; cell-number proportions",
    x = NULL, y = NULL
  ) +
  guides(colour = guide_legend(override.aes = list(shape = c(16, 15), size = 2.5))) +
  theme_fig3 +
  theme(
    axis.text.x = element_text(size = 8, colour = "#222222"),
    axis.text.y = element_text(size = 8, lineheight = 0.92, colour = "#222222"),
    panel.spacing.x = unit(2.5, "mm"),
    legend.position = "top",
    legend.margin = margin(t = 0, b = 0, unit = "pt"),
    legend.box.margin = margin(t = 0, b = 0, unit = "pt")
  )

# Panel c: positive expression in both references; preserve the no-extra-floor subset.
c_data <- crossref[crossref$min_top_mean == 0, , drop = FALSE]
c_data <- c_data[match(c(1, 1.2, 1.5, 2), c_data$min_dominance_ratio), , drop = FALSE]
if (nrow(c_data) != 4L || anyNA(c_data$min_dominance_ratio)) {
  stop("The min_top_mean == 0 source subset must contain all four locked ratio thresholds.")
}
c_data$threshold <- factor(
  c("All", "≥1.2", "≥1.5", "≥2.0"),
  levels = c("All", "≥1.2", "≥1.5", "≥2.0")
)
c_data$agreement_percent <- 100 * c_data$concordance
c_data$label <- sprintf("%.1f%%\n(n=%d)", c_data$agreement_percent, c_data$n_evaluable)
c_data$source_table <- basename(source_files[["cross_reference_concordance"]])
c_data$threshold_scope <- "min_top_mean == 0; positive expression in both references"

data_checks <- data.frame(
  check = c(
    "Panel A values match both frozen median-rho columns",
    "Panel B values match frozen CCC and MAE columns (MAE × 100 for percent display)",
    "Panel C values match the four min_top_mean == 0 rows",
    "Panel C locked count sequence is 220/341, 144/179, 103/111, 65/67",
    "External-cohort GSE108524 T-cell method median rho equals frozen source value"
  ),
  status = c(
    isTRUE(all.equal(a_data$rho, c(stability$method_effect_med, stability$reference_effect_med), tolerance = 0)),
    isTRUE(all.equal(b_data$value, c(calibration$CCC, 100 * calibration$MAE, 100 * rare$MAE_rare), tolerance = 0)),
    isTRUE(all.equal(c_data$concordance, c(0.645161290322581, 0.804469273743017, 0.927927927927928, 0.970149253731343), tolerance = 1e-14)),
    identical(c_data$n_concordant, c(220L, 144L, 103L, 65L)) && identical(c_data$n_evaluable, c(341L, 179L, 111L, 67L)),
    sum(stability$cohort == "GSE108524" & stability$cell_type == "T cell") == 1L &&
      abs(stability$method_effect_med[stability$cohort == "GSE108524" & stability$cell_type == "T cell"] - (-0.0580879434110543)) < 1e-14
  ),
  stringsAsFactors = FALSE
)
data_checks$status <- ifelse(data_checks$status, "PASS", "FAIL")
if (any(data_checks$status != "PASS")) stop("A plotted-data source check failed.")

p_c <- ggplot(c_data, aes(threshold, agreement_percent, group = 1)) +
  geom_line(colour = "#1B7837", linewidth = 0.7) +
  geom_point(shape = 21, size = 2.8, stroke = 0.7, fill = "white", colour = "#1B7837") +
  geom_text(
    aes(label = label), family = "Arial", size = 8 / .pt,
    lineheight = 0.95, vjust = -0.62, colour = "#222222"
  ) +
  scale_y_continuous(
    breaks = c(60, 70, 80, 90, 100), limits = c(55, 115),
    expand = expansion(mult = 0)
  ) +
  scale_x_discrete(expand = expansion(add = 0.22)) +
  labs(
    tag = "c",
    title = "Cross-reference dominant-compartment agreement",
    subtitle = "Positive expression in both references; no additional >0.05 expression floor",
    x = "Top-to-second expression ratio", y = "Agreement (%)"
  ) +
  theme_fig3 +
  theme(
    axis.text.x = element_text(size = 8, colour = "#222222"),
    axis.text.y = element_text(size = 8, colour = "#222222"),
    axis.title.x = element_text(size = 8, margin = margin(t = 1.5, unit = "mm")),
    axis.title.y = element_text(size = 8, margin = margin(r = 1.5, unit = "mm"))
  )

fig <- p_a / p_b / p_c + plot_layout(heights = c(0.95, 0.95, 0.92))
fig <- fig & theme(plot.background = element_rect(fill = "white", colour = NA))

width_mm <- 174
height_mm <- 198
base_name <- "Figure_3"
pdf_path <- file.path(output_dir, paste0(base_name, ".pdf"))
svg_path <- file.path(output_dir, paste0(base_name, ".svg"))
png_path <- file.path(output_dir, paste0(base_name, "_600dpi.png"))

font_match <- systemfonts::match_fonts("Arial")
if (nrow(font_match) < 1L || is.na(font_match$path[1]) || !nzchar(font_match$path[1])) {
  stop("Arial could not be resolved on this system.")
}

grDevices::cairo_pdf(
  filename = pdf_path, width = width_mm / 25.4, height = height_mm / 25.4,
  family = "Arial", bg = "white", onefile = TRUE
)
print(fig)
grDevices::dev.off()

svglite::svglite(
  file = svg_path, width = width_mm / 25.4, height = height_mm / 25.4,
  bg = "white", system_fonts = list(sans = "Arial")
)
print(fig)
grDevices::dev.off()

ragg::agg_png(
  filename = png_path, width = width_mm, height = height_mm,
  units = "mm", res = 600, background = "white"
)
print(fig)
grDevices::dev.off()

write.csv(a_data, file.path(plotdata_dir, "Figure_3_panel_a_plotted_data.csv"), row.names = FALSE)
write.csv(b_data, file.path(plotdata_dir, "Figure_3_panel_b_plotted_data.csv"), row.names = FALSE)
write.csv(c_data, file.path(plotdata_dir, "Figure_3_panel_c_plotted_data.csv"), row.names = FALSE)
write.csv(data_checks, file.path(qa_dir, "DATA_VALUE_CHECKS.csv"), row.names = FALSE)

source_manifest <- data.frame(
  source_table = basename(source_files),
  source_path = unname(source_files),
  output_copy = file.path("inputs", unname(source_files)),
  source_sha256 = vapply(unname(source_paths), function(path) {
    line <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE, stderr = TRUE)
    if (!length(line)) return(NA_character_)
    strsplit(trimws(line[1]), "[[:space:]]+")[[1]][1]
  }, character(1)),
  stringsAsFactors = FALSE
)
write.csv(source_manifest, file.path(source_manifest_dir, "figure_3_source_manifest.csv"), row.names = FALSE)

pdf_info <- system2("pdfinfo", shQuote(pdf_path), stdout = TRUE, stderr = TRUE)
page_line <- grep("^Pages:", pdf_info, value = TRUE)
size_line <- grep("^Page size:", pdf_info, value = TRUE)
if (length(page_line) != 1L || length(size_line) != 1L) stop("Could not inspect PDF page count or size.")
page_count <- as.integer(sub("^Pages:[[:space:]]*", "", page_line))
size_tokens <- regmatches(size_line, regexec("^Page size:[[:space:]]*([0-9.]+)[[:space:]]+x[[:space:]]+([0-9.]+)[[:space:]]+pts", size_line))[[1]]
if (length(size_tokens) != 3L) stop("Could not parse PDF page size.")
pdf_width_mm <- as.numeric(size_tokens[2]) * 25.4 / 72
pdf_height_mm <- as.numeric(size_tokens[3]) * 25.4 / 72
preview_prefix <- file.path(qa_dir, "Figure_3_pdf_render")
preview_status <- system2("pdftoppm", c("-f", "1", "-l", "1", "-r", "300", "-png", shQuote(pdf_path), shQuote(preview_prefix)))
if (!identical(preview_status, 0L)) stop("Could not render PDF QA preview.")
if (page_count != 1L) stop("PDF must have exactly one page.")
if (abs(pdf_width_mm - width_mm) > 0.15 || abs(pdf_height_mm - height_mm) > 0.15) {
  stop("PDF page dimensions do not match the requested physical size.")
}
pdf_text <- paste(system2("pdftotext", c("-layout", shQuote(pdf_path), "-"), stdout = TRUE, stderr = TRUE), collapse = "\n")
writeLines(pdf_text, file.path(qa_dir, "Figure_3_pdf_text.txt"), useBytes = TRUE)
required_tokens <- c(
  "Cross-configuration rank concordance", "GSE108524", "T cell", "0.06",
  "Method effect", "Reference effect", "Schwann", "tumor", "Myeloid", "Fibroblast",
  "Endothelial", "Pericyte", "CCC", "Overall MAE", "Rare-proportion MAE",
  "Donor-held-out synthetic-mixture recovery", "within-dataset holdout",
  "13 donors", "100 mixtures/donor", "64.5%", "80.4%", "92.8%", "97.0%",
  "341", "179", "111", "67"
)
missing_tokens <- required_tokens[!vapply(required_tokens, grepl, logical(1), x = pdf_text, fixed = TRUE)]
if (length(missing_tokens)) stop("PDF text check failed for: ", paste(missing_tokens, collapse = ", "))
if (grepl(input_root, pdf_text, fixed = TRUE)) stop("Private package input path found in PDF text.")
svg_text <- paste(readLines(svg_path, warn = FALSE), collapse = "\n")
if (!grepl('font-family: "Arial"', svg_text, fixed = TRUE)) stop("SVG text does not declare Arial.")
if (grepl(input_root, svg_text, fixed = TRUE)) stop("Private package input path found in SVG source.")

qa_summary <- data.frame(
  check = c(
    "PDF page count", "PDF width (mm)", "PDF height (mm)",
    "Panel A plotted rows", "Panel B plotted rows", "Panel C plotted rows",
    "GSE108524 T-cell method median rho", "Panel C source subset min_top_mean"
  ),
  result = c(
    page_count, sprintf("%.3f", pdf_width_mm), sprintf("%.3f", pdf_height_mm),
    nrow(a_data), nrow(b_data), nrow(c_data),
    format(stability$method_effect_med[stability$cohort == "GSE108524" & stability$cell_type == "T cell"], digits = 12),
    paste(unique(c_data$min_top_mean), collapse = ",")
  ),
  status = c("PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS"),
  stringsAsFactors = FALSE
)
write.csv(qa_summary, file.path(qa_dir, "BUILD_QA.csv"), row.names = FALSE)

message("Created PDF: ", pdf_path)
message("Created SVG: ", svg_path)
message("Created 600 dpi PNG: ", png_path)
message("Rendered PDF QA preview: ", paste0(preview_prefix, "-1.png"))
message("PDF dimensions: ", sprintf("%.1f × %.1f mm", pdf_width_mm, pdf_height_mm))
message("Arial resolution: ", font_match$path[1])
