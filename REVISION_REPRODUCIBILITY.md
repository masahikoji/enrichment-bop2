# Reviewer-revision reproducibility and release notes

Snapshot assembled 2026-10-25 from previously supplied R sources and CSV outputs. All binary-endpoint probabilities below use exact finite-state enumeration, not random simulation. Complex endpoints use Monte Carlo.

## What is included

- Unmodified base main and complex-endpoint pipeline structure, reorganized earlier to run relative to the GitHub root (`R/01...04`, `run_main.R`, `run_supplement.R`, `run_all.R`).
- R1-M1/M3: exact recalibration for additional `(p0,p1)` settings (R1 benchmark code).
- R1-M2: exact evaluation of reversed subgroup response ordering with fixed boundaries.
- R2-C2/C3: exact subgroup-only and claim-specific analyses with fixed boundaries. The code also evaluates an earlier exploratory theta-positive=0.30 scenario; manuscript Table S9 retains only theta-positive=0.40.
- R2-C4: standalone component calibration and exact global error-rate check for the adjusted BOP2 comparator; source uses the *largest feasible attained maximum error* criterion rather than minimizing distance to the proposed design.
- Main binary runtime and case-study calibration source code, plus their previously supplied output snapshots.

## How to reproduce

From the repository root, run one mode at a time:

```bash
Rscript run_main.R
Rscript run_supplement.R
Rscript run_revision.R r1-benchmarks
Rscript run_revision.R r1-reverse
Rscript run_revision.R r2-claims
Rscript run_revision.R r2-adjusted
Rscript run_revision.R main-timing
Rscript run_revision.R case-study
```

Set `BOP2_N_CORES=4` to limit workers when applicable. S4 Monte Carlo is slow and should use the prespecified full replication counts, not QUICK mode, for a manuscript reproduction. All new revision-run outputs go to **new directories** under `results/revision/` and will not overwrite the committed reference snapshots. For direct execution of R2 scripts, run from `R/revision/` so `source("R2_exact_common.R")` resolves.

R1 benchmarks default to the proposed-design gamma grid `seq(0.05,8.00,by=0.05)` used in the main manuscript. Supplied R-produced `ALL_summary.csv` and `ALL_selected_designs.csv` for all three benchmark combinations are committed under `reference_results/R1_benchmarks/primary_gamma8/` and designated as the reproducibility reference for Supplementary Table S7. An earlier **exploratory** extension to gamma upper bound 12 is retained separately under `reference_results/R1_benchmarks/exploratory_gamma12/` and can be repeated by setting `BOP2_PROPOSED_GAMMA_MAX=12`. Its selected positive gamma is 10.30 or 11.35 for the two added settings, compared with 7.80 or 7.95 for the prespecified 8-bound search; the reported 2-decimal operating characteristics are essentially unchanged. Do not mix parameter estimates or file versions from these distinct searches.

## Script / manuscript correspondence

| Item | Script | Snapshot |
|---|---|---|
| Published BOP2 binary Table I verification | `R/validation/01_reproduce_bop2_boundaries.R` | `reference_results/BOP2_validation/` |
| Main binary OC | `R/01_main_binary_exact.R` | generated under `results/main/` |
| Main tables | `R/02_main_tables.R` | generated under `tables/main/` |
| S4 complex endpoints | `R/03_supp_complex_type1.R` | generated under `results/supplement/` |
| S4 tables | `R/04_supp_tables.R` | generated under `tables/supplement/` |
| S5 sensitivity | `R/revision/R1M1_additional_settings_exact.R` | `reference_results/R1_benchmarks/primary_gamma8/ALL_summary.csv` and `ALL_selected_designs.csv` |
| S6 reversed ordering | `R/revision/R1M2_reverse_order_exact.R` | `reference_results/R1_reversed/` |
| S7 subgroup only | `R/revision/R2_Comments2_3_subgroup_only_exact.R` | `reference_results/R2_claims/` |
| S8 adjusted BOP2 | `R/revision/R2_Comment4_error_matched_BOP2.R` | `reference_results/R2_adjusted/` |
| S3 main calibration time | `R/revision/main_binary_boundary_calibration_timing.R` | `reference_results/timing/` |
| Prospective case-study calibration | `R/revision/case_study_boundary_calibration_timing.R` | `reference_results/case_study/` |

## Selected current numerical reference values

- Proposed point-global-null maximum across pi=0.4,0.5,0.6,0.7,0.8: 0.09894542782786911.
- Independent-comparator point-global-null maximum: 0.12633818863934865.
- Adjusted-comparator tau: 0.0874601406339 (reported rounded as approximately 0.08746); attained maximum 0.0989178353869672.
- Main 30-scenario mean PRN-any: proposed 83.45%, independent 78.76%, adjusted 77.33%. Proposed exceeds adjusted in 29 of 30 scenarios.
- Subgroup-only setting pi=0.4, theta_positive=0.4, theta_negative=0.05, theta_all=0.19: proposed PRN-positive 28.72%, PRN-all 4.63%.

These values are cross-checked against the committed reference CSVs using `python3 scripts/verify_reference_results.py`. The validator additionally checks the submitted gamma-8 benchmark results, the separate gamma-12 exploratory results, and the published BOP2 boundary reproduction snapshot. The independent benchmark CSVs were provided by the authors, not regenerated in the current environment.

## Methodological scope and unresolved verification

1. **Standalone BOP2 selection rule -- verified against the published numerical example:** Zhou, Lee and Yuan (2017, *Statistics in Medicine*, DOI 10.1002/sim.7338), Section 4, explicitly choose an attainable type I error closest to, but not above, the nominal level for discrete boundaries; the numerical study then considers power. The supplied implementation uses this rule: select the maximum attainable type I error <= cap, then maximize power among ties. The independently computed boundaries for the original paper's binary trial example are **(1,2,4,5,7,9,10)** at looks **(10,15,20,25,30,35,40)**, exactly the published Table I. Exact enumeration gives a point-null rejection probability 0.0974284191403399 and power 0.881487130225181, close to the publication's simulated 9.6% and 88.3%. This does **not** mean the selected design globally maximizes power among every candidate with type I error <= cap. We checked the published paper's described and tabulated numerical rule; we did not audit the authors' original web app code.
2. **S5 gamma grid -- primary/exploratory outputs distinguished:** the author-supplied R results for the **primary, manuscript-stated gamma upper bound 8** are now committed; exploratory upper-bound-12 outputs are separately labeled and must not be substituted for the primary parameter estimates. If the manuscript still describes all numerical specifications as unchanged, ensure it uses the gamma-8 primary results rather than the exploratory gamma-12 parameter estimates.
3. **S8 criterion:** only selection of the shared nominal level tau was changed in the R2-C4 code: among feasible candidate levels, maximize the *attained* point-global-null maximum. Standalone component ranking and all other numerical functions were left unchanged. The candidate-search CSV previously supplied shows this selects the same tau as the earlier comparison-with-proposed criterion. The R code with this adjustment has not been executed in the current environment (Rscript unavailable).
4. **Scope of error control:** all calibrations and reported maxima here refer to the *prespecified point global null* at the specified prevalence grid. This does not establish uniform strong FWER control over the full composite-null parameter space.
5. **No original patient data:** all files in `reference_results/` are generated analysis results. The dovitinib illustration is trial-inspired, not a reanalysis of participant-level records.

## Source provenance

- Existing GitHub tree was inspected (default branch `main`), and this package was based on the previously prepared repository-relative revision of that pipeline.
- New revision analyses were copied from supplied user files. `R2_exact_common.R` was the corrected 2026-10-24 version (entry_n=17, b_20=4).
- Code changes in this package are limited to repo-relative wrappers, output directory support, documentation, and the S8 tau selection criterion described above.
- A local Python integrity check was run against reference CSVs. A full R execution/Monte Carlo reproduction was not possible in this environment. Original BOP2 paper text and Table I were directly compared with the supplied stand-alone boundary reproduction and computed results. The user-supplied paper PDF is not included for copyright reasons.
