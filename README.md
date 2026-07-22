# Globally calibrated enrichment BOP2: reproducibility code

This repository reproduces the binary-endpoint analyses in the main manuscript and the complex-endpoint analyses in the Supplementary Material.

The implementation uses the **no-bridging rule**: when the all-comer futility boundary is crossed with fewer than 10 biomarker-positive patients available, the trial stops rather than accruing additional patients solely to reach the minimum subgroup sample size.

## Designs compared

### Proposed design

The all-comer and biomarker-positive boundary systems are jointly calibrated for the complete branching adaptive procedure. The calibration criterion controls

\[
\Pr(R_A \cup R_+)
\]

under the prespecified point global null over the biomarker-prevalence grid.

### Comparator

The comparator consists of independently calibrated conventional BOP2 components:

1. A conventional all-comer BOP2 is calibrated for analyses at 20, 30, and 40 patients.
2. For each attainable biomarker-positive entry sample size \(m=10,\ldots,30\), a conventional positive-population BOP2 is calibrated for the schedule consisting of \(m\), any subsequent scheduled analyses among 20 and 30 patients, and the final analysis at 40 patients.
3. For each standalone BOP2, candidates with type I error no greater than 0.10 are considered first; the attainable type I error closest to 0.10 is selected, and power is used to resolve ties.
4. The fixed standalone boundaries are then embedded in the no-bridging adaptive enrichment procedure. No constraint is imposed on the union probability \(\Pr(R_A\cup R_+)\) for the comparator.

## Requirements

- R 4.2 or later
- Base R packages for the analyses
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

Parallel workers can be specified as follows:

```bash
BOP2_N_CORES=4 Rscript run_all.R
```

## Main manuscript: binary endpoint

`R/01_main_binary_exact.R` performs exact recursive enumeration and writes results to `results/main/`.

Important outputs include:

- `exact_operating_characteristics.csv`
- `exact_type1_error.csv`
- `exact_power.csv`
- `selected_designs.csv`
- `proposed_boundary_table_allcomer.csv`
- `proposed_boundary_table_positive.csv`
- `comparator_allcomer_bop2_boundary.csv`
- `comparator_positive_bop2_library.csv`
- `enrichment_bop2_exact_results.xlsx` when an Excel-writing package is installed

`R/02_main_tables.R` converts these results to LaTeX tables in `tables/main/`.

## Supplementary Material: complex categorical endpoints

`R/03_supp_complex_type1.R` performs simulation-based calibration and independent validation for:

- nested efficacy;
- co-primary efficacy; and
- joint efficacy--toxicity monitoring.

The proposed design is globally calibrated. The comparator uses an independently calibrated conventional all-comer BOP2 and an entry-specific library of conventional positive-population BOP2 designs. The script also writes both proposed and comparator boundary outputs to `results/supplement/`.

`R/04_supp_tables.R` creates the corresponding LaTeX tables in `tables/supplement/`, including the standalone comparator boundary tables.

## Quick code checks

The reduced settings below are intended only to check code execution. They do not reproduce manuscript results.

```bash
BOP2_QUICK=1 Rscript run_main.R
BOP2_COMPLEX_QUICK=1 Rscript run_supplement.R
```

The complex-endpoint Monte Carlo replicate counts can be overridden:

```bash
BOP2_COMPLEX_N_CALIB=50000 BOP2_COMPLEX_N_VALID=100000 Rscript run_supplement.R
```

## Output mapping

See `MANUSCRIPT_OUTPUTS.md` for the mapping between manuscript tables and generated files. Generated result and table files are excluded by `.gitignore`; the directory placeholders are retained.
