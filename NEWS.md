## maicplus 0.1.2.9000

- Add vignette about proportional hazards
- Add continuous outcomes support for MAIC analysis
  - New effect measures: Mean Difference (MD), Standardized Mean Difference (SMD), Ratio of Means (RoM)
  - SMD methods: Cohen's d, Hedges' g, Glass's delta
  - New helper function `get_pseudo_ipd_continuous()` to generate pseudo IPD from aggregate continuous data
  - New helper function `lm_makeup()` to summarize linear model results
  - Support for both anchored and unanchored MAIC with continuous endpoints
  - Example datasets: `adlb_sat` and `adlb_twt` for continuous outcomes
  - Comprehensive vignettes for anchored and unanchored continuous analysis

## maicplus 0.1.2

- Changed the link function to binomial for RD and RR in binary analysis
- Updated vignette to recommend against bootstrap in unanchored analysis with citation
- Fix events percentage to properly incorporate weights in the report table for binary analysis
- Replace quantile function to use matrixStats package instead of desctools package
- Added testing PH assumptions vignette
- Fixed the issues with merging weights using match function

## maicplus 0.1.1

- bug fixes

## maicplus 0.1.0

- first CRAN release
