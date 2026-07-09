#!/usr/bin/env Rscript
# 20_prepare_submission_assets.R
# Prepare BMC-style supplementary/source-data assets from already generated result tables.
# This script does not recompute analyses and does not modify raw data.

options(stringsAsFactors = FALSE, warn = 1)

root <- Sys.getenv("PROJ_ROOT", getwd())
source_dir <- file.path(root, "results", "source_data")
supp_dir <- file.path(root, "results", "supplementary_tables")
dir.create(source_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)

copy_asset <- function(from, to_dir, to_name = basename(from)) {
  src <- file.path(root, from)
  if (!file.exists(src)) stop("missing source table: ", from)
  dst <- file.path(to_dir, to_name)
  ok <- file.copy(src, dst, overwrite = TRUE)
  if (!ok) stop("failed to copy: ", src, " -> ", dst)
  data.frame(source = from, output = sub(paste0("^", root, "/?"), "", dst))
}

manifest <- do.call(rbind, list(
  copy_asset("results/deg/deg_full_gene.csv", source_dir, "Fig1_DEG_full_gene.csv"),
  copy_asset("results/deg/deg_sig_fdr05_fc1.csv", source_dir, "Fig1_DEG_sig_fdr05_fc1.csv"),
  copy_asset("results/enrich/GO_BP_up.csv", source_dir, "Fig1_GO_BP_up.csv"),
  copy_asset("results/enrich/GO_BP_down.csv", source_dir, "Fig1_GO_BP_down.csv"),
  copy_asset("results/validation/hub_logFC_concordance.csv", source_dir, "Fig2_hub_logFC_concordance.csv"),
  copy_asset("results/scrna/celltype_composition.csv", source_dir, "Fig3_celltype_composition.csv"),
  copy_asset("results/scrna/gene_celltype_expression.csv", source_dir, "Fig3_gene_celltype_expression.csv"),
  copy_asset("results/pseudobulk/dominant_celltype_direction_concordance.csv", source_dir, "Fig4_direction_concordance.csv"),
  copy_asset("results/deconv/cell_fractions_GSE39645_music.csv", source_dir, "Fig4_GSE39645_MuSiC_fractions.csv"),
  copy_asset("results/deconv/deg_composition_adjusted.csv", source_dir, "Fig4_composition_adjusted_DEG.csv"),
  copy_asset("results/pseudobulk/intrinsic_vs_composition.csv", source_dir, "Fig4_intrinsic_vs_composition.csv"),
  copy_asset("results/gse216783/compartment_composition.csv", source_dir, "Fig5_GSE216783_compartment_composition.csv"),
  copy_asset("results/gse216783/candidate_compartment_localization.csv", source_dir, "Fig5_candidate_compartment_localization.csv"),
  copy_asset("results/gse216783/reference_comparison_deconv.csv", source_dir, "Fig5_reference_comparison_deconv.csv"),
  copy_asset("results/wgcna/module_trait_cor.csv", source_dir, "FigS1_module_trait_cor.csv"),
  copy_asset("results/hub/module_gene_metrics.csv", source_dir, "FigS1_module_gene_metrics.csv"),
  copy_asset("results/immune/ssgsea_scores.csv", source_dir, "FigS2_ssgsea_scores.csv"),
  copy_asset("results/immune/immune_diff_VS_vs_control.csv", source_dir, "FigS2_immune_diff_VS_vs_control.csv"),
  copy_asset("results/immune/signature_immune_cor.csv", source_dir, "FigS2_candidate_immune_cor.csv"),
  copy_asset("results/pseudobulk/pseudobulk_bulk_logFC_correlation.csv", source_dir, "FigS4_bulk_pseudobulk_correlation.csv"),
  copy_asset("results/pseudobulk/pseudobulk_bulk_logFC_sensitivity.csv", source_dir, "FigS4_bulk_pseudobulk_sensitivity.csv"),
  copy_asset("results/deconv/cell_fractions_GSE108524_music.csv", source_dir, "FigS5_GSE108524_MuSiC_fractions.csv"),
  copy_asset("results/deconv/cell_fraction_group_tests.csv", source_dir, "FigS5_cell_fraction_group_tests.csv"),
  copy_asset("results/deconv/method_agreement.csv", source_dir, "FigS5_method_agreement.csv"),
  copy_asset("results/wgcna/module_preservation.csv", source_dir, "FigS6_module_preservation.csv"),
  copy_asset("results/compartment_cross_evidence_map.csv", supp_dir, "TableS1_compartment_cross_evidence_map.csv")
))

write.csv(manifest, file.path(source_dir, "source_data_manifest.csv"), row.names = FALSE)
cat("Prepared", nrow(manifest), "source/supplementary assets\n")
cat("Source data:", source_dir, "\n")
cat("Supplementary tables:", supp_dir, "\n")
