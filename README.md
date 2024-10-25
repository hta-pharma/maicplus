<!-- markdownlint-disable MD013 MD033 -->

# maicplus <a href="https://hta-pharma.github.io/maicplus/"><img src="man/figures/logo.png" align="right" height="139" alt="maicplus website" /></a>

<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN
status](https://www.r-pkg.org/badges/version-last-release/maicplus)](https://www.r-pkg.org/badges/version-last-release/maicplus)
[![Code
Coverage](https://raw.githubusercontent.com/hta-pharma/maicplus/_xml_coverage_reports/data/main/badge.svg)](https://hta-pharma.github.io/maicplus/main/coverage-report/)

<!-- badges: end -->
<!-- markdownlint-enable MD013 MD033 -->

## Overview

The maicplus package provides tools to perform matching-adjusted indirect
comparison (MAIC) analysis. MAIC is a method used to adjust for differences
in baseline characteristics between treatment groups in indirect comparisons,
typically where patient-level data are available for one treatment but only
aggregate-level data are available for the comparator treatment.

## Features

- Perform MAIC analysis for time-to-event endpoints (e.g. overall survival) or binary endpoints (e.g. objective tumor response).

- Perform unanchored analysis when comparing single-arm trials or anchored analysis when there exists a common comparator between trials.

- User-friendly functions to facilitate the analysis process.

- Comprehensive documentation and examples to guide users through the analysis.

## Installation

You can install the latest version of `maicplus` on CRAN with:

```r
install.packages('maicplus')
```

or you can install the development version with:

```r
remotes::install_github("hta-pharma/maicplus")
```

## Tutorial

To learn how to use the `maicplus` R package, refer to the
[package website (hta-pharma.github.io/maicplus/)](https://hta-pharma.github.io/maicplus/).

## Bibliography

Chen G, Seo M, Antoniou M, Belleli R, Kalyvas C, Gravestock I. SA83 {Maicplus}:
An R Package to Support Analysis and Reporting of Matching Adjusted Indirect Treatment
Comparisons (MAIC) for HTA Dossiers. Value in Health. 2023;26(12):S558. [doi:10.1016/j.jval.2023.09.2992](https://doi.org/10.1016/j.jval.2023.09.2992)

Phillippo D, Ades T, Dias S, Palmer S, Abrams KR, Welton N.
NICE DSU Technical Support Document 18:
Methods for Population-Adjusted Indirect Comparisons in Submissions to NICE.
Vol 18. NICE Decision Support Unit; 2016.
