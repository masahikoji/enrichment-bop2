# Public release checklist (reviewer revision)

- [x] GitHub base structure retained, with extra sources under `R/revision/`.
- [x] Reviewer 1 and Reviewer 2 supplementary analysis code added.
- [x] Selected exact analysis CSVs committed separately under `reference_results/`.
- [x] Code's adjusted-comparator tau selection criterion no longer depends on the proposed design's realized error rate.
- [x] No participant-level clinical data were added.
- [x] Verify BOP2 comparator selection order against the original paper (Section 4 numerical-study convention) and reproduce its published binary Table I stopping boundaries. Original web-app internals were not audited.
- [x] Include author-supplied R1-M1 `ALL_summary.csv` and `ALL_selected_designs.csv` for gamma upper bound 8 (primary) and 12 (exploratory), stored separately.
- [ ] Make sure the final S5 text and tables use the gamma-8 primary results, consistent with the stated main search grid.
- [ ] Run the package's complete analyses in R, compare to reference snapshots, and preserve session information.
- [ ] Ensure manuscript, supplementary tables, reviewer letter, README, and scripts use identical calibration methods, labels, and values.
- [ ] Add a research software license appropriate to the authors/institutional requirements (not included because permission has not been specified).
- [ ] After these steps, tag a versioned public release and optionally attach a Zenodo DOI.

Do not represent this bundle as a fully validated manuscript release until the unchecked items have been resolved.
