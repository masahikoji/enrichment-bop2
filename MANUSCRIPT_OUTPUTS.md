# Manuscript output mapping

## Main manuscript

| Manuscript content | Generated file | Source script |
|---|---|---|
| Proposed boundary table | `tables/main/table_proposed_boundaries_compact.tex` | `R/02_main_tables.R` |
| Exact global type I error | `tables/main/table_type1_exact.tex` | `R/02_main_tables.R` |
| Power, \(\theta_+=0.4\) | `tables/main/table_power_comparison_exact_theta_pos_0p4.tex` | `R/02_main_tables.R` |
| Power, \(\theta_+=0.5\) | `tables/main/table_power_comparison_exact_theta_pos_0p5.tex` | `R/02_main_tables.R` |
| Power, \(\theta_+=0.6\) | `tables/main/table_power_comparison_exact_theta_pos_0p6.tex` | `R/02_main_tables.R` |
| Main-table input list | `tables/main/main_tables_input.tex` | `R/02_main_tables.R` |

## Additional binary-endpoint outputs

| Content | Generated file |
|---|---|
| Selected tuning parameters | `tables/main/table_selected_designs.tex` |
| Proposed operating characteristics | `tables/main/table_proposed_oc_exact_theta_pos_*.tex` |
| No-bridging enrichment-path results | `tables/main/table_proposed_enrichment_path_exact_theta_pos_*.tex` |
| Conventional all-comer comparator | `tables/main/table_comparator_allcomer_bop2.tex` |
| Entry-specific positive comparator library | `tables/main/table_comparator_positive_bop2_library.tex` |
| Supplemental binary-table input list | `tables/main/supplemental_tables_input.tex` |

## Supplementary Material: complex endpoints

| Manuscript content | Generated file | Source script |
|---|---|---|
| Selected complex-endpoint designs | `tables/supplement/table_complex_endpoint_selected_designs.tex` | `R/04_supp_tables.R` |
| Proposed complex-endpoint boundaries | `tables/supplement/table_complex_endpoint_boundaries_compact.tex` | `R/04_supp_tables.R` |
| Complex-endpoint global type I error | `tables/supplement/table_complex_endpoint_type1.tex` | `R/04_supp_tables.R` |
| Standalone all-comer comparator boundaries | `tables/supplement/table_complex_endpoint_standalone_allcomer_bop2.tex` | `R/04_supp_tables.R` |
| Positive comparator library: nested endpoint | `tables/supplement/table_complex_endpoint_positive_bop2_library_nested.tex` | `R/04_supp_tables.R` |
| Positive comparator library: co-primary endpoint | `tables/supplement/table_complex_endpoint_positive_bop2_library_coprimary.tex` | `R/04_supp_tables.R` |
| Positive comparator library: efficacy--toxicity endpoint | `tables/supplement/table_complex_endpoint_positive_bop2_library_efficacy_toxicity.tex` | `R/04_supp_tables.R` |
| Complex-table input list | `tables/supplement/supplemental_complex_endpoint_tables_input.tex` | `R/04_supp_tables.R` |

## Primary result files

Binary results are written to `results/main/`. Complex-endpoint results are written to `results/supplement/`. In particular, the complex-endpoint simulation writes:

- `complex_endpoint_proposed_boundaries.csv`;
- `complex_endpoint_comparator_boundaries.csv`;
- `complex_endpoint_standalone_allcomer_bop2.csv`; and
- `complex_endpoint_standalone_positive_bop2_library.csv`.

## Reviewer-requested analyses (S5--S8) and case-study outputs

The script names and committed reference output files are mapped in `REVISION_REPRODUCIBILITY.md`.

| Revised manuscript item | Code or committed reference result |
|---|---|
| S5 / Table S7: sensitivity to `(p0,p1)` | `R/revision/R1M1_additional_settings_exact.R`; `reference_results/R1_benchmarks/primary_gamma8/ALL_summary.csv` and `ALL_selected_designs.csv` (separate exploratory `gamma12/` outputs retained) |
| S6 / Table S8: reversed response ordering | `R/revision/R1M2_reverse_order_exact.R`; `reference_results/R1_reversed/R1M2_reverse_table_percent.csv` |
| S7 / Table S9: biomarker-positive-only efficacy | `R/revision/R2_Comments2_3_subgroup_only_exact.R`; `reference_results/R2_claims/R2_C2_subgroup_only_percent.csv` |
| S8 / Table S10: adjusted comparator, point-global-null error | `R/revision/R2_Comment4_error_matched_BOP2.R`; `reference_results/R2_adjusted/R2_C4_Point_global_null_percent.csv` |
| S8 / Table S11: 30 alternative settings | same script; `reference_results/R2_adjusted/R2_C4_Original_30_alternatives_percent.csv` |
| Case study: prospectively selected boundaries | `R/revision/case_study_boundary_calibration_timing.R`; `reference_results/case_study/` |
| Main runtime / S3 | `R/revision/main_binary_boundary_calibration_timing.R`; `reference_results/timing/` |

The R1-M1 `ALL_summary.csv` is not included because the actual R-run output for the chosen gamma grid was not available. The `reference_results/R1_benchmarks/` JSON files are *independent Python validation*, not R outputs.
