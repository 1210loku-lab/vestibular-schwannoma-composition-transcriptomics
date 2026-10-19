#!/usr/bin/env Rscript

# Rebuild the current JMN Fig. 1 / Fig. 2 panels as isolated, vector-first
# derivatives from the approved frozen figure-source tables.
# All source data are read-only; this script writes only to its output directory.

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
qa_dir <- file.path(output_root, "outputs", "figure_qa", "figures_1_2")
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
  dominant_compartment_concordance = "results/pseudobulk/dominant_celltype_direction_concordance.csv",
  fraction_group_tests = "results/deconv/cell_fraction_group_tests.csv",
  composition_adjusted_deg = "results/deconv/deg_composition_adjusted.csv",
  intrinsic_vs_composition = "results/pseudobulk/intrinsic_vs_composition.csv",
  hub_logfc_concordance = "results/validation/hub_logFC_concordance.csv",
  module_preservation = "results/wgcna/module_preservation.csv"
)
source_paths <- file.path(input_root, source_files)
missing_sources <- source_paths[!file.exists(source_paths)]
if (length(missing_sources)) stop("Missing frozen input table(s): ", paste(missing_sources, collapse = "; "))
tables <- lapply(source_paths, read.csv, check.names = FALSE)
names(tables) <- names(source_files)

theme_jmn <- theme_bw(base_family = "Arial", base_size = 8) +
  theme(
    text = element_text(family = "Arial", colour = "black", size = 8),
    axis.text = element_text(family = "Arial", colour = "black", size = 8),
    axis.title = element_text(family = "Arial", colour = "black", size = 8),
    plot.title = element_text(
      family = "Arial", colour = "black", size = 9, face = "bold", hjust = 0,
      margin = margin(l = 13, b = 1, unit = "pt")
    ),
    plot.subtitle = element_text(family = "Arial", colour = "grey25", size = 8, hjust = 0),
    plot.title.position = "plot",
    plot.tag = element_text(family = "Arial", colour = "black", size = 12, face = "bold"),
    plot.tag.position = c(0, 1),
    legend.title = element_text(family = "Arial", colour = "black", size = 8),
    legend.text = element_text(family = "Arial", colour = "black", size = 8),
    strip.text = element_text(family = "Arial", colour = "black", size = 8, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.25, colour = "grey90"),
    plot.margin = margin(t = 2, r = 3, b = 1, l = 3, unit = "mm"),
    legend.key.height = unit(3, "mm"),
    legend.key.width = unit(6, "mm")
  )

panel_tag_theme <- theme(plot.tag = element_text(family = "Arial", face = "bold", size = 12))
theme_set(theme_jmn)
tagger <- function(plot) plot + plot_annotation(tag_levels = "a", theme = panel_tag_theme)
write_plot_data <- function(x, name) {
  write.csv(x, file.path(plotdata_dir, name), row.names = FALSE, na = "")
}

# Figure 1, panel a: exact source-table concordance proportions and denominators.
dc <- tables$dominant_compartment_concordance
dca <- dc[dc$scope == "by_celltype" & dc$filter_type == "all", , drop = FALSE]
dca$cell_type <- ifelse(
  dca$dominant_celltype == "Schwann/肿瘤细胞", "Schwann/tumor",
  ifelse(dca$dominant_celltype == "巨噬/髓系", "Myeloid",
         ifelse(dca$dominant_celltype == "T细胞", "T cell",
                ifelse(dca$dominant_celltype == "成纤维", "Fibroblast",
                       ifelse(dca$dominant_celltype == "内皮", "Endothelial",
                              ifelse(dca$dominant_celltype == "周细胞", "Pericyte", dca$dominant_celltype)))))
)
dca$cell_type <- factor(dca$cell_type, levels = dca$cell_type[order(dca$direction_concordance)])
overall_row <- dc[dc$scope == "overall" & dc$filter_type == "all", , drop = FALSE]
if (nrow(dca) != 6L || nrow(overall_row) != 1L) stop("Fig. 1a source rows do not match the six-compartment/overall layout.")
overall_label <- sprintf("Overall %.1f%% (n=%d)", 100 * overall_row$direction_concordance, overall_row$n_direction)
p1a <- ggplot(dca, aes(direction_concordance, cell_type)) +
  geom_vline(xintercept = overall_row$direction_concordance, linetype = 2, colour = "grey35") +
  geom_col(fill = "#0072B2", width = 0.65) +
  geom_text(
    aes(label = sprintf("%.1f%% (n=%d)", 100 * direction_concordance, n_direction)),
    hjust = -0.05, family = "Arial", size = 8 / .pt
  ) +
  annotate("text", x = 1.45, y = 6.55, hjust = 1.02,
           family = "Arial", size = 8 / .pt, label = overall_label) +
  scale_x_continuous(limits = c(0, 1.5), breaks = c(0, 0.3, 0.6, 0.9),
                     labels = function(x) paste0(round(100 * x), "%")) +
  scale_y_discrete(expand = expansion(add = c(0.3, 1.0))) +
  labs(title = "Direction agreement by dominant compartment", x = "Direction agreement", y = NULL) +
  theme_jmn

# Figure 1, panel b: preserve all four source configurations and the current wording.
gt <- tables$fraction_group_tests
cell_order <- c("Schwann/tumor", "Myeloid", "T cell", "Fibroblast", "Endothelial", "Pericyte")
zh2en <- c(
  "Schwann/肿瘤细胞" = "Schwann/tumor", "巨噬/髓系" = "Myeloid", "T细胞" = "T cell",
  "成纤维" = "Fibroblast", "内皮" = "Endothelial", "周细胞" = "Pericyte"
)
std_cell <- function(x) {
  z <- unname(zh2en[as.character(x)])
  ifelse(is.na(z), as.character(x), z)
}
stars <- function(p) {
  ifelse(is.na(p), "NE", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns"))))
}
gt$cell_type <- factor(std_cell(gt$cell_type), levels = rev(cell_order))
gt$method <- factor(gt$method, levels = c("MuSiC", "NNLS"))
gt$cohort <- factor(gt$cohort, levels = c("GSE39645", "GSE108524"))
gt$configuration <- factor(
  paste(gt$cohort, gt$method, sep = "\n"),
  levels = c("GSE39645\nMuSiC", "GSE108524\nMuSiC", "GSE39645\nNNLS", "GSE108524\nNNLS")
)
gt$mark <- stars(gt$padj)
p1b <- ggplot(gt, aes(configuration, cell_type, fill = median_shift)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(aes(label = mark), family = "Arial", size = 8 / .pt, fontface = "bold") +
  scale_fill_gradient2(
    low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0,
    breaks = c(-0.2, -0.1, 0, 0.1, 0.2),
    name = "Median VS-minus-control score",
    guide = guide_colourbar(direction = "horizontal", title.position = "top", title.hjust = 0.5,
                            barwidth = unit(32, "mm"), barheight = unit(2, "mm"))
  ) +
  labs(title = "Estimated-composition group contrasts", x = NULL, y = NULL) +
  theme_jmn +
  theme(
    axis.text.x = element_text(size = 8, angle = 60, hjust = 1, vjust = 1, lineheight = 0.93),
    axis.text.y = element_text(size = 8),
    legend.position = "top", legend.justification = c(1, 0.5),
    legend.box.margin = margin(t = 0, b = 0, unit = "pt"),
    panel.grid = element_blank()
  )

# Figure 1, panel c: display every row in the frozen adjusted-contrast table.
ca <- tables$composition_adjusted_deg
ca$display <- ifelse(
  ca$classification == "attenuated", "Lost one or both thresholds",
  ifelse(ca$classification == "retained", "Threshold support retained", "other genes")
)
ca$display <- factor(ca$display, levels = c("Lost one or both thresholds", "Threshold support retained", "other genes"))
ca_cols <- c("Lost one or both thresholds" = "#56B4E9", "Threshold support retained" = "#D55E00", "other genes" = "#BDBDBD")
p1c <- ggplot(ca, aes(original_logFC, adjusted_logFC, colour = display)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey45") +
  geom_point(size = 0.35, alpha = 0.45) +
  scale_colour_manual(values = ca_cols, name = NULL, breaks = names(ca_cols)) +
  labs(
    title = "Composition-adjusted differential expression",
    subtitle = "326/1,410 original DEGs retained both thresholds",
    x = expression(Original~log[2]~FC), y = expression(Adjusted~log[2]~FC)
  ) +
  theme_jmn + theme(legend.position = "bottom", legend.box = "horizontal")

# Figure 1, panel d: category names are mapped to the current author-edited wording.
ivc <- tables$intrinsic_vs_composition
class_map <- c(
  "composition-consistent" = "Composition-aligned",
  both = "Mixed-evidence",
  "cell-intrinsic (candidate)" = "Directionally supported\nwithin-compartment",
  "ambiguous/undetermined" = "Unresolved"
)
ivc$display <- unname(class_map[ivc$class])
if (anyNA(ivc$display)) stop("Unmapped evidence category in the frozen Fig. 1d source table.")
class_order <- c("Unresolved", "Directionally supported\nwithin-compartment", "Mixed-evidence", "Composition-aligned")
cc <- as.data.frame(table(factor(ivc$display, levels = class_order)), stringsAsFactors = FALSE)
names(cc) <- c("class", "count")
cc$class <- factor(as.character(cc$class), levels = class_order)
cc$percent <- 100 * cc$count / sum(cc$count)
class_cols <- c(
  "Composition-aligned" = "#0072B2", "Mixed-evidence" = "#7B3294",
  "Directionally supported\nwithin-compartment" = "#D55E00", "Unresolved" = "#BDBDBD"
)
p1d <- ggplot(cc, aes(count, class, fill = class)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%d (%.1f%%)", count, percent)), hjust = -0.05,
            family = "Arial", size = 8 / .pt) +
  scale_fill_manual(values = class_cols, guide = "none", breaks = names(class_cols)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.45))) +
  labs(title = "Exploratory screening categories of bulk DEGs", x = "Number of discovery DEGs", y = NULL) +
  theme_jmn + theme(axis.text.y = element_text(size = 8))

fig1 <- tagger((p1a | p1b) / (p1c | p1d) + plot_layout(heights = c(1, 1)))

# Figure 2, panel a: all 346 source rows, retaining the author-displayed summary.
val <- tables$hub_logfc_concordance
val$module <- factor(val$module, levels = c("turquoise", "yellow"))
module_colors <- c(turquoise = "#009E73", yellow = "#E69F00", other = "#999999")
if (nrow(val) != 346L) stop("Fig. 2a source row count no longer matches the author-displayed fixed set (n=346).")
p2a <- ggplot(val, aes(discovery_logFC, validation_logFC, colour = module)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.4, colour = "grey45") +
  geom_point(size = 1.0, alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.65, colour = "black") +
  scale_colour_manual(values = module_colors[c("turquoise", "yellow")], name = "Module") +
  annotate("label", x = -Inf, y = Inf, hjust = -0.04, vjust = 1.1,
           family = "Arial", size = 8 / .pt, linewidth = 0.2,
           label = "Pearson r = 0.883\nDirection agreement = 93.4%") +
  labs(title = "Cross-platform fold-change reproducibility",
       x = expression(Discovery~log[2]~FC), y = expression(GSE108524~log[2]~FC)) +
  coord_cartesian(clip = "off") + theme_jmn + theme(legend.position = "bottom")

# Figure 2, panel b: retain the source script's top-14 selection and current display order.
topg <- head(val[order(-abs(val$discovery_logFC)), , drop = FALSE], 14)
top_order <- as.character(topg$gene[order(topg$discovery_logFC)])
topg$gene <- factor(topg$gene, levels = rev(top_order))
p2b <- ggplot(topg, aes(y = gene)) +
  geom_segment(aes(x = discovery_logFC, xend = validation_logFC, yend = gene),
               colour = "grey65", linewidth = 0.5) +
  geom_point(aes(x = discovery_logFC, shape = "Discovery"), size = 2.0, colour = "#7B3294") +
  geom_point(aes(x = validation_logFC, shape = "External"), size = 2.0, colour = "#009E73") +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  scale_shape_manual(values = c(Discovery = 16, External = 17), name = NULL) +
  labs(title = "Largest discovery effects across cohorts", x = expression(log[2]~fold~change), y = NULL) +
  theme_jmn +
  theme(axis.text.y = element_text(family = "Arial", face = "italic", size = 8), legend.position = "bottom")

# Figure 2, panel c: same module-preservation rows, order, thresholds, and colors.
mp <- tables$module_preservation
mp$module_class <- ifelse(mp$module %in% c("turquoise", "yellow"), mp$module, "other")
mp$module <- factor(mp$module, levels = mp$module[order(mp$Zsummary)])
p2c <- ggplot(mp, aes(Zsummary, module, colour = module_class)) +
  geom_vline(xintercept = c(2, 10), linetype = 2, colour = c("#0072B2", "#009E73"), linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = Zsummary, yend = module), colour = "grey75", linewidth = 0.5) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = module_colors, guide = "none") +
  labs(title = "Cross-platform preservation of discovery modules",
       subtitle = "Dashed lines: Zsummary = 2 and 10", x = "Preservation Zsummary", y = "Module") +
  theme_jmn

fig2 <- tagger((p2a | p2b) / p2c + plot_layout(heights = c(1, 0.9)))

# Source-exact plotted data and transparent provenance copies.
write_plot_data(dca, "Figure_1_panel_a_plotted_data.csv")
write_plot_data(gt, "Figure_1_panel_b_plotted_data.csv")
write_plot_data(ca, "Figure_1_panel_c_plotted_data.csv")
write_plot_data(cc, "Figure_1_panel_d_plotted_data.csv")
write_plot_data(val, "Figure_2_panel_a_plotted_data.csv")
write_plot_data(topg, "Figure_2_panel_b_plotted_data.csv")
write_plot_data(mp, "Figure_2_panel_c_plotted_data.csv")

source_manifest <- data.frame(
  source_table = basename(source_files),
  source_path = unname(source_files),
  output_copy = file.path("inputs", unname(source_files)),
  source_sha256 = vapply(source_paths, function(path) {
    line <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE, stderr = TRUE)
    if (!length(line)) return(NA_character_)
    strsplit(trimws(line[1]), "[[:space:]]+")[[1]][1]
  }, character(1)),
  stringsAsFactors = FALSE
)
write.csv(source_manifest, file.path(source_manifest_dir, "figures_1_2_source_manifest.csv"), row.names = FALSE)

data_checks <- data.frame(
  check = c(
    "Fig. 1a percentages and denominators match six by-celltype/all source rows",
    "Fig. 1b includes the 24 source rows (6 cell types × 4 cohort-method configurations)",
    "Fig. 1c retains every frozen adjusted-DE table row and current classification mapping",
    "Fig. 1d category counts match 688/378/52/292",
    "Fig. 2a uses all 346 fixed-set source rows",
    "Fig. 2b top-14 selected gene set matches current author labels",
    "Fig. 2c includes all ten frozen module-preservation rows"
  ),
  status = c(
    nrow(dca) == 6L && nrow(overall_row) == 1L && setequal(as.integer(dca$n_direction), c(154L, 131L, 235L, 69L, 319L, 454L)),
    nrow(gt) == 24L,
    nrow(ca) == 17328L && all(ca$display %in% names(ca_cols)),
    identical(as.integer(cc$count[match(c("Composition-aligned", "Mixed-evidence", "Directionally supported\nwithin-compartment", "Unresolved"), cc$class)]), c(688L, 378L, 52L, 292L)),
    nrow(val) == 346L,
    identical(sort(as.character(topg$gene)), sort(c("FABP4", "MFAP5", "PRRX1", "CCDC80", "PLA2G2A", "TMOD2", "LOC653653", "PDZD2", "CAB39L", "LGI4", "ARHGEF26", "RASSF4", "NCAM2", "SOX2"))),
    nrow(mp) == 10L
  ),
  stringsAsFactors = FALSE
)
data_checks$status <- ifelse(data_checks$status, "PASS", "FAIL")
if (any(data_checks$status != "PASS")) stop("A frozen-source or plotted-data check failed.")
write.csv(data_checks, file.path(qa_dir, "DATA_VALUE_CHECKS.csv"), row.names = FALSE)

font_match <- systemfonts::match_fonts("Arial")
if (nrow(font_match) < 1L || is.na(font_match$path[1]) || !nzchar(font_match$path[1])) stop("Arial could not be resolved.")

export_one <- function(plot, stem, height_mm) {
  width_mm <- 174
  pdf_path <- file.path(output_dir, paste0(stem, ".pdf"))
  svg_path <- file.path(output_dir, paste0(stem, ".svg"))
  png_path <- file.path(output_dir, paste0(stem, "_600dpi.png"))

  grDevices::cairo_pdf(
    filename = pdf_path, width = width_mm / 25.4, height = height_mm / 25.4,
    family = "Arial", bg = "white", onefile = TRUE
  )
  print(plot)
  grDevices::dev.off()

  svglite::svglite(
    file = svg_path, width = width_mm / 25.4, height = height_mm / 25.4,
    bg = "white", system_fonts = list(sans = "Arial")
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(filename = png_path, width = width_mm, height = height_mm,
                units = "mm", res = 600, background = "white")
  print(plot)
  grDevices::dev.off()

  info <- system2("pdfinfo", shQuote(pdf_path), stdout = TRUE, stderr = TRUE)
  page_line <- grep("^Pages:", info, value = TRUE)
  size_line <- grep("^Page size:", info, value = TRUE)
  if (length(page_line) != 1L || length(size_line) != 1L) stop("Could not inspect PDF dimensions for ", stem)
  page_count <- as.integer(sub("^Pages:[[:space:]]*", "", page_line))
  size_tokens <- regmatches(size_line, regexec("^Page size:[[:space:]]*([0-9.]+)[[:space:]]+x[[:space:]]+([0-9.]+)[[:space:]]+pts", size_line))[[1]]
  if (length(size_tokens) != 3L) stop("Could not parse PDF page size for ", stem)
  pdf_w <- as.numeric(size_tokens[2]) * 25.4 / 72
  pdf_h <- as.numeric(size_tokens[3]) * 25.4 / 72
  if (page_count != 1L || abs(pdf_w - width_mm) > 0.3 || abs(pdf_h - height_mm) > 0.3) {
    stop("Unexpected PDF page count or physical dimensions for ", stem)
  }
  preview_prefix <- file.path(qa_dir, paste0(stem, "_pdf_render"))
  preview_status <- system2("pdftoppm", c("-f", "1", "-l", "1", "-r", "300", "-png", shQuote(pdf_path), shQuote(preview_prefix)))
  if (!identical(preview_status, 0L)) stop("Could not render PDF QA preview for ", stem)
  preview <- paste0(preview_prefix, "-1.png")
  pdf_text <- paste(system2("pdftotext", c("-layout", shQuote(pdf_path), "-"), stdout = TRUE, stderr = TRUE), collapse = "\n")
  writeLines(pdf_text, file.path(qa_dir, paste0(stem, "_pdf_text.txt")), useBytes = TRUE)
  if (grepl(input_root, pdf_text, fixed = TRUE)) stop("Private package input path found in PDF text for ", stem)
  svg_text <- paste(readLines(svg_path, warn = FALSE), collapse = "\n")
  if (!grepl('font-family: "Arial"', svg_text, fixed = TRUE)) stop("SVG does not declare Arial for ", stem)
  if (grepl(input_root, svg_text, fixed = TRUE)) stop("Private package input path found in SVG source for ", stem)
  list(
    stem = stem, pdf = pdf_path, svg = svg_path, png = png_path,
    preview = preview, width_mm = pdf_w, height_mm = pdf_h,
    pdf_text = pdf_text
  )
}

fig1_out <- export_one(fig1, "Figure_1", 175)
fig2_out <- export_one(fig2, "Figure_2", 160)

required <- list(
  Fig1 = c(
    "Direction agreement by dominant compartment", "Overall 69.2% (n=1362)",
    "Estimated-composition group contrasts", "Median VS-minus-control score",
    "Composition-adjusted differential expression", "326/1,410 original DEGs retained both thresholds",
    "Lost one or both thresholds", "Threshold support retained",
    "Exploratory screening categories of bulk DEGs", "Composition-aligned", "Mixed-evidence",
    "Directionally supported", "within-compartment", "Unresolved", "688 (48.8%)", "378 (26.8%)",
    "52 (3.7%)", "292 (20.7%)", "Schwann/tumor"
  ),
  Fig2 = c(
    "Cross-platform fold-change reproducibility", "Pearson r = 0.883", "Direction agreement = 93.4%",
    "GSE108524", "FABP4", "MFAP5", "PRRX1", "CCDC80", "PLA2G2A", "TMOD2", "LOC653653",
    "PDZD2", "CAB39L", "LGI4", "ARHGEF26", "RASSF4", "NCAM2", "SOX2", "External",
    "Cross-platform preservation of discovery modules", "Dashed lines: Zsummary = 2 and 10"
  )
)
for (out in list(fig1_out, fig2_out)) {
  panel <- if (identical(out$stem, "Figure_1")) "Fig1" else "Fig2"
  missing <- required[[panel]][!vapply(required[[panel]], grepl, logical(1), x = out$pdf_text, fixed = TRUE)]
  if (length(missing)) stop(panel, " PDF text check failed: ", paste(missing, collapse = ", "))
}

qa_summary <- data.frame(
  figure = c("Fig1", "Fig2"),
  pdf_file = c(basename(fig1_out$pdf), basename(fig2_out$pdf)),
  page_width_mm = sprintf("%.3f", c(fig1_out$width_mm, fig2_out$width_mm)),
  page_height_mm = sprintf("%.3f", c(fig1_out$height_mm, fig2_out$height_mm)),
  png_resolution_dpi = 600L,
  font_family = "Arial",
  visual_review = "R-rendered PDF preview created; inspect before integration",
  status = "PASS",
  stringsAsFactors = FALSE
)
write.csv(qa_summary, file.path(qa_dir, "BUILD_QA.csv"), row.names = FALSE)

message("Created Fig. 1: ", fig1_out$pdf)
message("Created Fig. 2: ", fig2_out$pdf)
message("SVG and 600-dpi PNG companions created for both figures.")
message("PDF dimensions (mm): Fig. 1 ", sprintf("%.1f × %.1f", fig1_out$width_mm, fig1_out$height_mm),
        "; Fig. 2 ", sprintf("%.1f × %.1f", fig2_out$width_mm, fig2_out$height_mm))
message("Arial resolved to: ", font_match$path[1])
