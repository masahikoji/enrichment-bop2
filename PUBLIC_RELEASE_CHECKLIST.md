# Public release checklist

- [x] No local file-system paths are embedded in the code.
- [x] No names, email addresses, affiliations, or other personal identifiers are included.
- [x] Main analysis outputs are generated under `results/main/` and corresponding LaTeX tables under `tables/main/`.
- [x] Complex-endpoint outputs are generated under `results/supplement/` and corresponding LaTeX tables under `tables/supplement/`.
- [x] Main calibration uses the single prespecified working alternative and does **not** use an alpha-lower selection band.
- [x] The three detailed binary enrichment-path tables in the Supplementary Material are generated and mapped in `MANUSCRIPT_OUTPUTS.md`.
- [x] Quick-check outputs are isolated from full-analysis outputs.
- [ ] Add the chosen software license before public release.
- [ ] Run the full analyses from a clean checkout and compare all generated manuscript values with the submitted manuscript.
- [ ] Record the R version and platform from the generated `sessionInfo.txt` files.
