# Globally calibrated enrichment BOP2: reproducibility code

This repository reproduces the numerical results and LaTeX tables reported in the main manuscript and supplementary material.

See `MANUSCRIPT_OUTPUTS.md` for a direct mapping between manuscript table labels and generated files.

## Requirements

- R 4.2 or later
- Base R packages only for the analyses
- Optional: `writexl` or `openxlsx` for consolidated Excel workbooks

Run commands from the repository root.

```bash
Rscript run_main.R
Rscript run_supplement.R
```

To reproduce everything:

```bash
Rscript run_all.R
```

Parallel workers can be set before execution:

```bash
BOP2_N_CORES=4 Rscript run_all.R
```

## Main manuscript

`R/01_main_binary_exact.R` performs exact recursive enumeration for the binary endpoint and writes results to `results/main/`.

`R/02_main_tables.R` creates:

- `tables/main/table_proposed_boundaries_compact.tex` (`tab:proposed_boundaries_compact`)
- `tables/main/table_type1_exact.tex` (`tab:type1_exact`)
- `tables/main/table_power_comparison_exact_theta_pos_0p4.tex`
- `tables/main/table_power_comparison_exact_theta_pos_0p5.tex`
- `tables/main/table_power_comparison_exact_theta_pos_0p6.tex`
- `tables/main/main_tables_input.tex`

## Supplementary material

`R/03_supp_complex_type1.R` performs simulation-based calibration and independent type I error evaluation for nested efficacy, co-primary efficacy, and joint efficacy-toxicity endpoints. Results are written to `results/supplement/`.

`R/04_supp_tables.R` creates:

- `tables/supplement/table_complex_endpoint_boundaries_compact.tex` (`tab:complex_endpoint_boundaries_compact`)
- `tables/supplement/table_complex_endpoint_type1.tex` (`tab:complex_endpoint_type1`)
- `tables/supplement/supplement_tables_input.tex`

## Quick checks

The reduced settings below are intended only for code checks and do not reproduce manuscript results.

```bash
BOP2_QUICK=1 Rscript run_main.R
BOP2_COMPLEX_QUICK=1 Rscript run_supplement.R
```

Generated result and table files are ignored by Git by default. Remove the corresponding entries from `.gitignore` if selected outputs should be version-controlled.
