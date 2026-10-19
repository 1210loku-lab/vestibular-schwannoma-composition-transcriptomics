#!/usr/bin/env Rscript

# Publication-size derivatives from frozen precomputed analysis tables only; no model fitting.
ROOT <- Sys.getenv("VS_OUTPUT_ROOT", unset = ".")
ANALYSIS <- file.path(ROOT, "outputs")
FIGURES <- file.path(ROOT, "figures")
dir.create(FIGURES, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

# Typography contract: Arial throughout; lowercase panel tags are 12 pt bold,
# panel titles are 9 pt bold, subtitles are 8 pt regular. No letter-spacing or
# tracking override is set; glyph spacing uses the font's default metrics.
theme_publication <- function(base_size = 10) {
  theme_bw(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(family = "Arial", size = 9, colour = "#222222"),
      axis.title = element_text(family = "Arial", size = 9),
      axis.text = element_text(family = "Arial", size = 9),
      strip.text = element_text(family = "Arial", face = "bold", size = 9),
      plot.title = element_text(family = "Arial", face = "bold", size = 9),
      plot.subtitle = element_text(family = "Arial", face = "plain", size = 8),
      plot.caption = element_blank(),
      plot.tag = element_text(family = "Arial", face = "bold", size = 12),
      legend.title = element_text(family = "Arial", size = 9),
      legend.text = element_text(family = "Arial", size = 9),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.grid.major = element_line(colour = "grey90", linewidth = .25),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "#666666", fill = NA, linewidth = .25),
      axis.line = element_blank(),
      axis.ticks = element_line(colour = "#666666", linewidth = .25),
      legend.position = "bottom",
      plot.margin = margin(6, 7, 6, 7)
    )
}

save_formats <- function(plot, stem, width, height) {
  svglite(paste0(stem, ".svg"), width = width, height = height)
  print(plot)
  dev.off()
  cairo_pdf(paste0(stem, ".pdf"), width = width, height = height, family = "Arial")
  print(plot)
  dev.off()
  agg_png(paste0(stem, ".png"), width = width, height = height, units = "in", res = 600)
  print(plot)
  dev.off()
}

read_output_csv <- function(file) read.csv(file.path(ANALYSIS, file), check.names = FALSE,
                                        stringsAsFactors = FALSE)

# S15: fresh design and effect sensitivity, arranged for the 6.5-in final page width.
eff <- read_output_csv("Figure_S15_effect_threshold_source.csv")
eff <- eff[eff$model_class == "composition" & eff$analysis_set == "full_31v9", , drop = FALSE]
if (nrow(eff) != 12L) stop("S15 expects the 12 accepted full-cohort specifications")
eff$label <- factor(eff$label, levels = rev(eff$label))
ref_col <- c("MuSiC" = "#7B3294", "NNLS" = "#3B6FB6", "GSE216783 NNLS" = "#D55E00")
eff$reference_label <- factor(eff$reference_label, levels = names(ref_col))
eff$pc_label <- factor(paste0("PC1–", eff$pc_n), levels = paste0("PC1–", 2:5))
p15a <- ggplot(eff, aes(x = retained_percent, y = pc_label, colour = reference_label)) +
  geom_vline(xintercept = 100, colour = "#777777", linetype = "dashed", linewidth = .4) +
  geom_point(size = 2.4) +
  facet_wrap(~reference_label, nrow = 1, scales = "free_y") +
  scale_colour_manual(values = ref_col, drop = FALSE) +
  labs(title = "Joint-threshold retention across 12 model specifications", x = "Baseline DEGs retaining both thresholds (%)", y = "PC specification") +
  theme_publication() + theme(axis.text.y = element_text(size = 9), strip.text = element_text(size = 9)) +
  guides(colour = "none")

primary <- eff[eff$original_main_MuSiC_PC1_5, , drop = FALSE]
if (nrow(primary) != 1L) stop("S15 expects one primary MuSiC PC1-5 row")
state <- data.frame(
  state = factor(c("Retained both", "Lost FDR only", "Lost effect only", "Lost both"),
    levels = c("Retained both", "Lost FDR only", "Lost effect only", "Lost both")),
  count = c(primary$both_thresholds_retained, primary$only_FDR_threshold_lost,
            primary$only_effect_threshold_lost, primary$both_thresholds_lost)
)
state$label <- ifelse(state$count >= 20L, paste0(state$count, " (", sprintf("%.1f", 100 * state$count / primary$baseline_fresh_DEG_n), "%)"), "")
p15b <- ggplot(state, aes(x = "Primary MuSiC PC1-5", y = count, fill = state)) +
  geom_col(width = .62, colour = "white", linewidth = .4) +
  geom_text(aes(label = label), position = position_stack(vjust = .5), size = 3.2,
            family = "Arial", colour = "#202020") +
  coord_flip() +
  scale_fill_manual(values = c("Retained both" = "#D55E00", "Lost FDR only" = "#56B4E9",
    "Lost effect only" = "#7B3294", "Lost both" = "#BDBDBD")) +
  labs(title = "Threshold-state accounting in the fixed 1,410-gene set", subtitle = "119 direction reversals: 117 threshold-lost; 2 retained", x = NULL, y = "Baseline DEGs (n)", fill = NULL) +
  theme_publication() + theme(axis.text.y = element_text(size = 9), legend.position = "bottom",
    legend.title = element_blank(), legend.text = element_text(size = 9),
    legend.key.width = grid::unit(.26, "in"), legend.box = "vertical") +
  guides(fill = guide_legend(title = NULL, nrow = 2, byrow = TRUE))
if (sum(state$count) != 1410L) stop("S15 primary threshold-state counts must sum to 1,410")

se <- read_output_csv("Figure_S15_SE_ratio_source.csv")
se$pc_n <- as.integer(sub(".*PC1-", "", se$label))
se$pc_label <- factor(paste0("PC1–", se$pc_n), levels = paste0("PC1–", 2:5))
se$reference_label <- factor(se$reference_label, levels = names(ref_col))
se$precision <- factor(se$precision, levels = c("Raw SE ratio", "Moderated SE ratio"))
p15c <- ggplot(se, aes(x = estimate, y = pc_label, colour = reference_label, shape = precision)) +
  geom_vline(xintercept = 1, colour = "#777777", linetype = "dashed", linewidth = .4) +
  geom_point(position = position_dodge(width = .55), size = 2) +
  scale_colour_manual(values = ref_col, drop = FALSE) +
  facet_wrap(~reference_label, nrow = 1, scales = "free_y") +
  labs(title = "Raw and moderated standard-error ratios", x = "Median gene-wise SE ratio vs baseline", y = NULL, shape = "SE type") +
  theme_publication() + theme(axis.text.y = element_text(size = 9), legend.position = "bottom",
    legend.title = element_text(size = 9), legend.text = element_text(size = 9),
    legend.key.width = grid::unit(.25, "in")) +
  guides(colour = "none", shape = guide_legend(nrow = 1, byrow = TRUE))

gate <- read_output_csv("Figure_S15_gate_source.csv")
gate <- gate[gate$analysis_set == "full_31v9" | gate$analysis_set == "tissue_31v8", , drop = FALSE]
gate$reference_public <- ifelse(gate$analysis_set == "tissue_31v8", "Tissue-only", gate$reference_label)
gate$reference_public <- factor(gate$reference_public,
  levels = c(names(ref_col), "Tissue-only"))
gate$pc_label <- factor(paste0("PC1–", gate$pc_n), levels = paste0("PC1–", 2:5))
gate$pc_label[gate$analysis_set == "tissue_31v8"] <- "PC1–5"
gate$pc_label <- factor(gate$pc_label, levels = paste0("PC1–", 2:5))
gate_col <- c(ref_col, "Tissue-only" = "#777777")
p15d <- ggplot(gate, aes(x = group_VIF, y = pc_label, colour = reference_public,
                         shape = accepted_by_original_numeric_gate)) +
  geom_vline(xintercept = 10, colour = "#555555", linetype = "dashed", linewidth = .45) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = gate_col, drop = FALSE) +
  facet_wrap(~reference_public, nrow = 1, scales = "free_y",
    labeller = labeller(reference_public = c("GSE216783 NNLS" = "GSE216783\nNNLS"))) +
  scale_shape_manual(values = c("FALSE" = 4, "TRUE" = 16),
                     labels = c("FALSE" = "Gate warning", "TRUE" = "Gate passed")) +
  labs(title = "Group-VIF diagnostics", x = "Group VIF (prespecified limit 10)", y = "PC specification", shape = "Design gate") +
  theme_publication() + theme(axis.text.y = element_text(size = 9), legend.position = "bottom",
    legend.title = element_text(size = 9), legend.text = element_text(size = 9),
    legend.key.width = grid::unit(.25, "in")) +
  guides(colour = "none", shape = guide_legend(nrow = 1, byrow = TRUE))
p15 <- (p15a / p15b / p15c / p15d) + plot_layout(heights = c(1.2, .9, 1.2, 1.4)) +
  plot_annotation(tag_levels = "a") &
  theme_publication() & theme(plot.margin = margin(10, 12, 10, 12), legend.box = "vertical")
save_formats(p15, file.path(FIGURES, "Figure_S15"), 7.2, 10.7)

# S16: fixed 346-gene scatter and matched point/interval estimates.
sc <- read_output_csv("Figure_S16_scatter_source.csv")
sc$contrast_label <- factor(sc$contrast_label,
  levels = c("All VS (27)", "Sporadic (10)", "NF2-associated (17)"))
lims <- range(c(sc$discovery_logFC, sc$logFC), finite = TRUE)
p16a <- ggplot(sc, aes(x = discovery_logFC, y = logFC, colour = contrast_label)) +
  geom_abline(slope = 1, intercept = 0, colour = "#555555", linetype = "dashed", linewidth = .4) +
  geom_point(alpha = .5, size = .8, stroke = 0) +
  facet_wrap(~contrast_label, nrow = 1) +
scale_colour_manual(values = c("All VS (27)" = "#009E73", "Sporadic (10)" = "#0072B2", "NF2-associated (17)" = "#D55E00")) +
  coord_equal(xlim = lims, ylim = lims, expand = TRUE) +
  labs(title = "External effects for the fixed 346-gene set", x = "Discovery group+batch effect", y = "External group-only difference") +
  theme_publication() + theme(axis.text = element_text(size = 9), strip.text = element_text(size = 9)) +
  guides(colour = "none")

ci <- read_output_csv("Figure_S16_interval_source.csv")
ci <- ci[ci$discovery_model == "fresh_group_batch" & ci$metric %in% c("Direction match", "Pearson r", "Spearman ρ", "Lin CCC"), , drop = FALSE]
ci$comparison_label <- factor(ci$comparison_label,
  levels = c("All VS", "Sporadic", "NF2"), labels = c("All VS", "Sporadic", "NF2-related"))
ci$metric <- factor(ci$metric, levels = c("Direction match", "Pearson r", "Spearman ρ", "Lin CCC"))
if (nrow(sc) != 1038L || anyNA(sc$contrast_label) || any(table(sc$contrast_label) != 346L))
  stop("S16 scatter must contain 346 fixed genes in each of the three exact contrasts")
if (nrow(ci) != 12L || anyNA(ci$comparison_label) || any(table(ci$metric) != 3L) ||
    any(table(ci$comparison_label) != 4L))
  stop("S16 intervals must contain four metrics and one row per exact comparison in each metric")
p16b <- ggplot(ci, aes(x = comparison_label, y = point_estimate,
                       ymin = percentile_2.5, ymax = percentile_97.5,
                       colour = comparison_label)) +
  geom_hline(yintercept = 0, colour = "#777777", linetype = "dashed", linewidth = .35) +
  geom_linerange(linewidth = .65) +
  geom_point(size = 1.25) +
  facet_wrap(~metric, ncol = 2) +
  scale_colour_manual(values = c("All VS" = "#009E73", "Sporadic" = "#0072B2", "NF2-related" = "#D55E00")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Sample-resampled agreement metrics", subtitle = "500 fixed sample-stratified draws", x = NULL, y = "Estimate (95% percentile interval)") +
  theme_publication() + theme(axis.text.x = element_text(size = 9)) +
  guides(colour = "none")
p16 <- p16a / p16b + plot_layout(heights = c(1, .85)) +
  plot_annotation(tag_levels = "a") &
  theme_publication() & theme(plot.margin = margin(10, 12, 10, 12))
save_formats(p16, file.path(FIGURES, "Figure_S16"), 7.2, 8.4)

# S17: donor-level expression using only the saved English reader fields.
d <- read_output_csv("donor_Figure_S17_source.csv")
d <- d[, c("gene", "lineage", "donor_id", "group", "n_cells", "log2_normalized_count_plus_1"), drop = FALSE]
d$gene <- factor(d$gene, levels = c("NRG1", "L1CAM", "FCGBP", "PRRX1"))
d$group <- factor(d$group, levels = c("Nerve control", "Vestibular schwannoma (VS)"))
d$lineage_public <- ifelse(d$gene %in% c("NRG1", "L1CAM"), "Schwann/tumor-associated annotation",
  ifelse(d$gene == "FCGBP", "Myeloid", "Fibroblast"))
s17sum <- aggregate(log2_normalized_count_plus_1 ~ gene + lineage + group, d,
                    function(x) c(mean = mean(x), sd = sd(x)))
s17sum$mean <- s17sum$log2_normalized_count_plus_1[, "mean"]
s17sum$sd <- s17sum$log2_normalized_count_plus_1[, "sd"]
s17sum$lineage <- factor(s17sum$lineage, levels = unique(d$lineage))
gene_lineage <- c(NRG1 = "Combined Schwann/tumor-associated annotation",
                  L1CAM = "Combined Schwann/tumor-associated annotation",
                  FCGBP = "Myeloid annotation", PRRX1 = "Fibroblast annotation")
if (nrow(d) != 36L || any(table(d$gene, d$group) != matrix(rep(c(2L, 7L), 4L), nrow = 4L, byrow = TRUE)))
  stop("S17 requires four predefined genes with two nerve and seven VS donors each")
make_s17_gene_panel <- function(gene) {
  dd <- d[d$gene == gene, , drop = FALSE]
  if (length(unique(dd$lineage)) != 1L || nrow(dd) != 9L)
    stop(paste("S17 gene-lineage mapping is not one 9-donor pair:", gene))
  ss <- aggregate(log2_normalized_count_plus_1 ~ group, dd,
                  function(x) c(mean = mean(x), sd = sd(x)))
  ss$mean <- ss$log2_normalized_count_plus_1[, "mean"]
  ss$sd <- ss$log2_normalized_count_plus_1[, "sd"]
  gene_title <- parse(text = sprintf("bolditalic('%s')", gene))[[1]]
  ggplot(dd, aes(x = group, y = log2_normalized_count_plus_1, colour = group)) +
    geom_col(data = ss, aes(x = group, y = mean, fill = group, colour = group),
             inherit.aes = FALSE, width = .52, alpha = .45, linewidth = .55) +
    geom_errorbar(data = ss, aes(x = group, ymin = mean - sd, ymax = mean + sd, colour = group),
                  inherit.aes = FALSE, width = .15, linewidth = .55) +
    geom_point(position = position_jitter(width = .075, height = 0, seed = 20261015),
               size = 2.0, alpha = .95) +
    scale_x_discrete(labels = c("Nerve control" = "Nerve\nn=2", "Vestibular schwannoma (VS)" = "VS\nn=7")) +
    scale_colour_manual(values = c("Nerve control" = "#0072B2", "Vestibular schwannoma (VS)" = "#D55E00")) +
    scale_fill_manual(values = c("Nerve control" = "#0072B2", "Vestibular schwannoma (VS)" = "#D55E00")) +
    labs(title = gene_title, subtitle = gene_lineage[[gene]], x = NULL,
         y = "log2(DESeq2-normalized count + 1)") +
    theme_publication() + theme(legend.position = "none", plot.title = element_text(family = "Arial", face = "bold.italic", size = 9),
      plot.subtitle = element_text(family = "Arial", face = "plain", size = 8), axis.text.x = element_text(size = 9),
      axis.text.y = element_text(size = 9)) + guides(colour = "none", fill = "none")
}
pp17 <- lapply(c("NRG1", "L1CAM", "FCGBP", "PRRX1"), make_s17_gene_panel)
p17 <- (pp17[[1]] + pp17[[2]] + pp17[[3]] + pp17[[4]] + plot_layout(ncol = 2)) +
  plot_annotation(tag_levels = "a") &
  theme_publication() & theme(plot.margin = margin(8, 10, 8, 10))
save_formats(p17, file.path(FIGURES, "Figure_S17"), 7.2, 7.4)

# S18: frozen repeat-level simulation outputs; no additional model fits.
r_base <- read_output_csv("simulation_repeat_gene_class_metrics.csv")
r_strong <- read_output_csv("simulation_strong_repeat_gene_class_metrics.csv")
r <- rbind(r_base, r_strong)
r$scenario[r$scenario == "separation"] <- "moderate"
r <- r[r$scenario %in% c("overlap", "moderate", "strong"), , drop = FALSE]
if (nrow(r) != 1080L || anyNA(r$scenario) || anyNA(r$route) || anyNA(r$gene_class))
  stop("S18 requires all 1,080 rows across three scenarios, three routes and six gene classes")
if (any(table(r$scenario, r$route, r$gene_class) != 20L) || any(table(r$scenario) != 360L))
  stop("S18 requires 20 repeat-level observations per scenario-route-gene-class cell")
routes <- c("naive", "oracle_true_composition", "estimated_music_composition")
route_labels <- c("Naive", "Oracle", "Estimated\nMuSiC")
names(route_labels) <- routes
scenario_labels <- c("overlap" = "Overlap", "moderate" = "Moderate", "strong" = "Strong")
route_cols <- c("Naive" = "#999999", "Oracle" = "#009E73", "Estimated\nMuSiC" = "#7B3294")
route_fills <- c("Naive" = "#999999", "Oracle" = "#009E73", "Estimated\nMuSiC" = "#7B3294")
make_s18_panel <- function(dat, cls, value_col, title, ylab, percent = FALSE) {
  dat <- dat[dat$gene_class == cls & dat$route %in% routes, , drop = FALSE]
  if (nrow(dat) != 180L || any(table(dat$scenario, dat$route) != 20L))
    stop(paste("S18 panel does not contain 20 rows per scenario/route:", cls))
  dat$route_label <- factor(unname(route_labels[dat$route]), levels = unname(route_labels))
  dat$scenario_label <- factor(unname(scenario_labels[dat$scenario]), levels = scenario_labels)
  dat$value <- dat[[value_col]]
  if (percent) dat$value <- 100 * dat$value
  sm <- aggregate(value ~ scenario_label + route_label, dat,
                  function(x) c(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE)))
  sm$mean <- sm$value[, "mean"]
  sm$sd <- sm$value[, "sd"]
  sm$value <- NULL
  ggplot(dat, aes(x = route_label, y = value, colour = route_label, fill = route_label)) +
    geom_col(data = sm, aes(x = route_label, y = mean, fill = route_label, colour = route_label), inherit.aes = FALSE,
             width = .62, linewidth = .5) +
    geom_errorbar(data = sm, aes(x = route_label, ymin = mean - sd, ymax = mean + sd, colour = route_label),
                  inherit.aes = FALSE, width = .16, linewidth = .45) +
    geom_point(position = position_jitter(width = .1, height = 0, seed = 20261015),
               alpha = .42, size = 1.1) +
    facet_wrap(~scenario_label, nrow = 1) +
    scale_colour_manual(values = route_cols) + scale_fill_manual(values = route_fills) +
    labs(title = title, x = NULL, y = ylab) +
  theme_publication() + theme(axis.text.x = element_text(size = 9),
                       strip.text = element_text(size = 9), plot.title = element_text(family = "Arial", face = "bold", size = 9),
                       plot.subtitle = element_text(family = "Arial", face = "plain", size = 8)) +
    guides(colour = "none", fill = "none")
}
p18a <- make_s18_panel(r, "composition_only", "false_call_rate_beta_zero",
  "Conditional beta-zero calls", "Threshold-call rate (%)", TRUE)
p18b <- make_s18_panel(r, "both", "sensitivity_nonzero_beta",
  "Both-effect sensitivity", "Sensitivity (%)", TRUE)
p18c <- make_s18_panel(r, "expression_only", "median_raw_se",
  "Expression-only standard errors", "Median raw SE", FALSE)
p18 <- (p18a / p18b / p18c) + plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme_publication() & theme(plot.margin = margin(10, 12, 10, 12))
save_formats(p18, file.path(FIGURES, "Figure_S18"), 7.2, 8.7)

message("Final supplementary figures written from frozen CSV inputs only")
message("Inputs: Figure_S15_*.csv, Figure_S16_*.csv, donor_Figure_S17_source.csv, simulation_repeat_gene_class_metrics.csv, simulation_strong_repeat_gene_class_metrics.csv")
