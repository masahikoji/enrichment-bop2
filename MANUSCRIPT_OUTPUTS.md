# Manuscript output map

## Main manuscript

Run:

```bash
Rscript run_main.R
```

| Manuscript label | Generated file | Source analysis |
|---|---|---|
| `tab:proposed_boundaries_compact` | `tables/main/table_proposed_boundaries_compact.tex` | Exact binary-endpoint calibration |
| `tab:type1_exact` | `tables/main/table_type1_exact.tex` | Exact binary-endpoint type I error |
| `tab:power_comparison_exact_0p4` | `tables/main/table_power_comparison_exact_theta_pos_0p4.tex` | Exact binary-endpoint power, theta-positive = 0.4 |
| `tab:power_comparison_exact_0p5` | `tables/main/table_power_comparison_exact_theta_pos_0p5.tex` | Exact binary-endpoint power, theta-positive = 0.5 |
| `tab:power_comparison_exact_0p6` | `tables/main/table_power_comparison_exact_theta_pos_0p6.tex` | Exact binary-endpoint power, theta-positive = 0.6 |

The main numerical outputs are written to `results/main/`.

## Supplementary material: detailed binary enrichment path

These tables are produced by the same command, `Rscript run_main.R`.

| Supplementary label | Generated file | Source analysis |
|---|---|---|
| `tab:proposed_enrichment_path_exact_0p4` | `tables/main/table_proposed_enrichment_path_exact_theta_pos_0p4.tex` | Exact binary enrichment-path operating characteristics, theta-positive = 0.4 |
| `tab:proposed_enrichment_path_exact_0p5` | `tables/main/table_proposed_enrichment_path_exact_theta_pos_0p5.tex` | Exact binary enrichment-path operating characteristics, theta-positive = 0.5 |
| `tab:proposed_enrichment_path_exact_0p6` | `tables/main/table_proposed_enrichment_path_exact_theta_pos_0p6.tex` | Exact binary enrichment-path operating characteristics, theta-positive = 0.6 |

`tables/main/supplemental_tables_input.tex` inputs these three tables.

## Supplementary material: complex categorical endpoints

Run:

```bash
Rscript run_supplement.R
```

| Supplementary label | Generated file | Source analysis |
|---|---|---|
| `tab:complex_endpoint_boundaries_compact` | `tables/supplement/table_complex_endpoint_boundaries_compact.tex` | Complex-endpoint simulation calibration |
| `tab:complex_endpoint_type1` | `tables/supplement/table_complex_endpoint_type1.tex` | Independent complex-endpoint type I error simulation |

The supplementary numerical outputs are written to `results/supplement/`.
