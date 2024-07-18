## maicplus

`maicplus` is an open-source R package that facilitates performing Matching Adjusted Indirect Comparison (MAIC) analysis. This package is designed to handle endpoints of interest that are either time-to-event (e.g., overall survival) or binary (e.g., objective tumor response).

## Installation

You can install the development version of `maicplus` from GitHub with:

```R
# install.packages("devtools")
devtools::install_github("hta-pharma/maicplus")
```

## Overview
The maicplus package provides tools to perform MAIC analysis. MAIC is a method used to adjust for differences in baseline characteristics between treatment groups in indirect comparisons, typically where patient-level data are available for one treatment but only aggregate-level data are available for the comparator treatment. This package specifically supports endpoints that are time-to-event or binary.

## Features
* Perform MAIC analysis for time-to-event endpoints (e.g., overall survival).
* Perform MAIC analysis for binary endpoints (e.g., objective tumor response).
* User-friendly functions to facilitate the analysis process.
* Comprehensive documentation and examples to guide users through the analysis.

## Usage
* Time-to-Event Endpoint Example (use code for example)
* Binary Endpoint Example (use code for example)

## Documentation
Detailed documentation for each function is available within the package. You can access the documentation using the ? operator in R. Additionally, comprehensive vignettes are provided to help you get started and understand the application of the package: 
```R
vignette("maicplus").
```

## License
This project is licensed under the MIT License - see the (put the LICENSE here) file for details.

## Contact
If you have any questions or feedback, please get in touch with us at [email@example.com].
