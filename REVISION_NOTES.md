# Revision notes

This version replaces the earlier componentwise comparator implementation.

## Main binary analysis

- Uses the no-bridging enrichment rule.
- Calibrates the conventional all-comer BOP2 independently on the 20, 30, and 40 patient schedule.
- Calibrates a separate conventional biomarker-positive BOP2 for every attainable entry sample size from 10 through 30.
- Applies the conventional BOP2 selection order: attainable type I error closest to and not above 0.10, followed by maximum power.
- Embeds the fixed standalone comparator boundaries in the adaptive procedure and evaluates PRN-any afterward.
- Generates dedicated comparator boundary and calibration files.

## Complex-endpoint analysis

- Uses the same no-bridging rule.
- Retains global calibration for the proposed design.
- Replaces the previous componentwise joint-grid comparator with independently calibrated conventional BOP2 components.
- Generates proposed boundaries, comparator boundaries, the standalone all-comer calibration, and all entry-specific positive-population calibrations.

## Tables

- Labels the comparator as independent conventional BOP2 components.
- Adds standalone comparator boundary tables for the binary and complex-endpoint analyses.
