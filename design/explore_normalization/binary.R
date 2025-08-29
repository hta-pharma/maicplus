# setwd("~/GitHub/maicplus")
source("design/explore_normalization/normalize_weights.R")

library(maicplus)
library(dplyr)

data(centered_ipd_sat)
data(adrs_sat)

centered_colnames <- c("AGE", "AGE_SQUARED", "SEX_MALE", "ECOG0", "SMOKE", "N_PR_THER_MEDIAN")
centered_colnames <- paste0(centered_colnames, "_CENTERED")

# get dummy binary pseudo IPD
pseudo_adrs <- get_pseudo_ipd_binary(
  binary_agd = data.frame(
    ARM = "B",
    RESPONSE = c("YES", "NO"),
    COUNT = c(280, 120)
  ),
  format = "stacked"
)

weighted_data <- estimate_weights(
  data = centered_ipd_sat,
  centered_colnames = centered_colnames
)
weighted_data_SW1 <- weighted_data_SWN <- weighted_data_SWESS <- weighted_data
weighted_data_SW1$data$weights <- normalize_weights(weighted_data$data$weights, method = "SW1")
weighted_data_SWN$data$weights <- normalize_weights(weighted_data$data$weights, method = "SWN")
weighted_data_SWESS$data$weights <- normalize_weights(weighted_data$data$weights, method = "SWESS")


### OW
result_OW <- maic_unanchored(
  weights_object = weighted_data,
  ipd = adrs_sat,
  pseudo_ipd = pseudo_adrs,
  trt_ipd = "A",
  trt_agd = "B",
  normalize_weight = FALSE,
  endpoint_type = "binary",
  endpoint_name = "Binary Endpoint",
  eff_measure = "OR",
  # binary specific args
  binary_robust_cov_type = "HC3"
)

### SW1
result_SW1 <- maic_unanchored(
  weights_object = weighted_data_SW1,
  ipd = adrs_sat,
  pseudo_ipd = pseudo_adrs,
  trt_ipd = "A",
  trt_agd = "B",
  normalize_weight = FALSE,
  endpoint_type = "binary",
  endpoint_name = "Binary Endpoint",
  eff_measure = "OR",
  # binary specific args
  binary_robust_cov_type = "HC3"
)

### SWN
result_SWN <- maic_unanchored(
  weights_object = weighted_data_SWN,
  ipd = adrs_sat,
  pseudo_ipd = pseudo_adrs,
  trt_ipd = "A",
  trt_agd = "B",
  normalize_weight = FALSE,
  endpoint_type = "binary",
  endpoint_name = "Binary Endpoint",
  eff_measure = "OR",
  # binary specific args
  binary_robust_cov_type = "HC3"
)

### SWESS
result_SWESS <- maic_unanchored(
  weights_object = weighted_data_SWESS,
  ipd = adrs_sat,
  pseudo_ipd = pseudo_adrs,
  trt_ipd = "A",
  trt_agd = "B",
  normalize_weight = FALSE,
  endpoint_type = "binary",
  endpoint_name = "Binary Endpoint",
  eff_measure = "OR",
  # binary specific args
  binary_robust_cov_type = "HC3"
)

# results are exactly the same
result_OW$inferential$summary
result_SW1$inferential$summary
result_SWN$inferential$summary
result_SWESS$inferential$summary
