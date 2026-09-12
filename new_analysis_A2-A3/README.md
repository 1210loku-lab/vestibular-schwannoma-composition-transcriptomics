# External localization and deconvolution robustness analyses

This directory contains the analyses added for manuscript version 1.2.0.

`A2_redo.py` processes GSE250061 from the three GEO 10x HDF5 matrices, applies the stated quality-control and marker-based broad-lineage annotation workflow, and writes the per-transcript localization table and supporting plots. `A2_concordance_stats.py` calculates exact and top-two concordance, a bootstrap confidence interval, a label-permutation expectation and the confusion matrix. Run the scripts from this directory after placing the three GEO matrices under `gse/` and the transcript map as `candidate_priority_matrix.csv`.

`A3_part1_refxmethod.R` produces the two-reference by two-estimator robustness matrix and adds MuSiC with GSE230375 as a sensitivity configuration. `A3_part2_benchmark.R` performs the donor-held-out synthetic-mixture benchmark. `A3_part2b_calibration.R` writes mixture-level estimates and the CCC, MAE, calibration and rare-proportion summaries. These scripts expect the public-project inputs described in the repository README and use `PROJ_ROOT` when the project is not the current working directory.

`35_build_cbc_final_figure5.R` reconstructs the vector Figure 5 from `reference_by_method_stability.csv`, `known_mixture_calibration.csv`, `known_mixture_rare_lowprop.csv` and `global_compartment_concordance_summary.csv`. Scientific conclusions should be based on the CSV source tables rather than values read from the rendered figure.
