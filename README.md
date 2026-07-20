# Globally calibrated enrichment BOP2: reproducibility code

This repository reproduces the numerical results and LaTeX tables reported in the main manuscript and supplementary material for the globally calibrated adaptive-enrichment BOP2 design.

See `MANUSCRIPT_OUTPUTS.md` for the mapping between manuscript table labels and generated files.

## Requirements

- R 4.2 or later
- Base and recommended R packages only for the analyses
- Optional: `writexl` or `openxlsx` for consolidated Excel workbooks

Run all commands from the repository root.

```bash
Rscript run_main.R
Rscript run_supplement.R
```

To reproduce everything:

```bash
Rscript run_all.R
```

Parallel workers may be set before execution:

```bash
BOP2_N_CORES=4 Rscript run_all.R
```

## Main manuscript and binary-endpoint supplement

`R/01_main_binary_exact.R` performs exact recursive enumeration for the binary endpoint and writes results to `results/main/`.

The proposed design is selected from **all** candidate boundary pairs satisfying

```text
maximum global type I error over the prevalence grid <= 0.10
```

There is no lower type I error selection band. Among feasible candidates, the selection criteria are applied in this order:

1. largest average PRN-any under the single working alternative `(theta_positive, theta_negative) = (0.40, 0.20)` over the prevalence grid;
2. largest minimum PRN-any under that working alternative;
3. largest maximum global type I error, subject to remaining at or below 0.10.

The fixed selected design is then evaluated over `theta_positive` in `{0.40, 0.50, 0.60}` and `theta_negative` in `{0.20, 0.10}`.

`R/02_main_tables.R` creates the five main-manuscript tables and the three detailed binary enrichment-path tables in the Supplementary Material. Outputs are written to `tables/main/`.

## Complex categorical endpoints

`R/03_supp_complex_type1.R` performs simulation-based calibration and independent type I error evaluation for nested efficacy, co-primary efficacy, and joint efficacy-toxicity endpoints. Results are written to `results/supplement/`.

The proposed complex-endpoint design requires the maximum one-sided 95% Wilson upper confidence bound for PRN-any over the prevalence grid to be at most 0.10. Among feasible candidates, the least conservative design is selected by average, then minimum, then maximum estimated PRN-any. The componentwise comparator separately constrains the Wilson upper bounds for PRN-all and PRN-positive and selects the feasible pair closest to `(0.10, 0.10)` in squared Euclidean distance.

`R/04_supp_tables.R` creates the complex-endpoint boundary and type I error tables in `tables/supplement/`.

## Quick checks

The reduced settings below are intended only for code checks and do not reproduce manuscript results. Quick-check outputs are stored in `quick_test` subdirectories and do not overwrite full-analysis outputs.

```bash
BOP2_QUICK=1 Rscript run_main.R
BOP2_COMPLEX_QUICK=1 Rscript run_supplement.R
```

## Output policy

Generated result and table files are ignored by Git by default. Remove or modify the corresponding `.gitignore` entries when selected outputs should be version-controlled.

No local user paths or personal identifiers are embedded in the analysis scripts.
