# Globally calibrated enrichment BOP2: reproducibility code

This repository contains analysis code for the binary-endpoint main manuscript, complex categorical endpoints, and reviewer-requested revision analyses. The original main/supplementary entry points are retained. A separate reviewer-revision entry point and reference CSVs have been added.

**Revision update (2026-10-25):** the original BOP2 numerical-study selection rule was verified against Zhou, Lee and Yuan (2017, Section 4) and its binary stopping boundaries in Table I were reproduced exactly. User-supplied R outputs for both the prespecified gamma upper bound 8 and exploratory upper bound 12 are now included. The **gamma <= 8** results are the primary manuscript reference, consistent with the stated candidate grid. See [REVISION_REPRODUCIBILITY.md](REVISION_REPRODUCIBILITY.md). Scientific code was not re-run here; the remaining items are identified in the release checklist.

The implementation uses the **no-bridging rule**: when the all-comer futility boundary is crossed with fewer than 10 biomarker-positive patients available, the trial stops rather than accruing additional patients solely to reach the minimum subgroup sample size.

## Designs compared

### Proposed design

The all-comer and biomarker-positive boundary systems are jointly calibrated for the complete branching adaptive procedure. The calibration criterion controls

\[
\Pr(R_A \cup R_+)
\]

under the prespecified point global null over the biomarker-prevalence grid.

### Comparator

The comparator used in the current manuscript consists of independently calibrated BOP2 components. The *implemented* selection order follows the **numerical-study convention in Zhou et al. (2017, Section 4)**. Their Section 2.3 gives a broad power-optimality description; Section 4 additionally specifies choosing an attainable type I error rate closest to, but not greater than, the nominal value for discrete boundaries. The selection steps are:

1. A conventional all-comer BOP2 is calibrated for analyses at 20, 30, and 40 patients.
2. For each attainable biomarker-positive entry sample size \(m=10,\ldots,30\), a conventional positive-population BOP2 is calibrated for the schedule consisting of \(m\), any subsequent scheduled analyses among 20 and 30 patients, and the final analysis at 40 patients.
3. For each standalone BOP2, candidates with type I error no greater than 0.10 are considered first; the attainable type I error closest to 0.10 is selected, and power is used to resolve ties, as in the original numerical study. This lexicographic convention does **not** claim to maximize power over all candidate designs with type I error <= 0.10.
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

## Reviewer-requested revision analyses (R1/R2)

Reproducibility instructions, exact analysis scenarios, output mapping, and caveats:
**[REVISION_REPRODUCIBILITY.md](REVISION_REPRODUCIBILITY.md)**.

Run individual analyses from the repository root:

```bash
Rscript run_revision.R r1-benchmarks
Rscript run_revision.R r1-reverse
Rscript run_revision.R r2-claims
Rscript run_revision.R r2-adjusted
Rscript run_revision.R main-timing
Rscript run_revision.R case-study
```

For the explored extended gamma grid in the R1 sensitivity analysis:

```bash
BOP2_PROPOSED_GAMMA_MAX=12 Rscript run_revision.R r1-benchmarks
```

Precomputed results are stored in `reference_results/` as data-only snapshots. New output directories are created below `results/revision/` and ignored by Git. This repository does not include actual participant-level data.

Run the lightweight cross-file result checks without R:

```bash
python3 scripts/verify_reference_results.py
```

## Original BOP2 validation

To verify the BOP2 numerical-study selection rule against the published binary example and the current all-comer schedule:

```bash
Rscript R/validation/01_reproduce_bop2_boundaries.R
```

The script writes under `results/revision/bop2_validation/`; previously computed snapshots are committed in `reference_results/BOP2_validation/`. The published stopping boundaries are reproduced as **(1, 2, 4, 5, 7, 9, 10)** for interim/final looks (10, 15, 20, 25, 30, 35, 40). The cited paper is Zhou H, Lee JJ, Yuan Y (2017), *Statistics in Medicine*, DOI: [10.1002/sim.7338](https://doi.org/10.1002/sim.7338). The article PDF itself is not redistributed here.

## R1 benchmark reproducibility snapshots

`reference_results/R1_benchmarks/primary_gamma8/` contains the author-supplied R outputs for the main (0.20, 0.40) and two additional benchmark combinations, using the stated proposed-design gamma grid 0.05 to 8.00. The results exactly reproduce the manuscript-level reported rounded values. `reference_results/R1_benchmarks/exploratory_gamma12/` contains separately labeled exploratory extended-grid outputs; **these are not the designated Table S7 reference values**. To reproduce the primary benchmark outputs, run `Rscript run_revision.R r1-benchmarks` without setting `BOP2_PROPOSED_GAMMA_MAX`.
