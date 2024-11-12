#
# unanchored example using maic_unanchored for binary outcome
#

data(centered_ipd_sat)
data(adrs_sat)

centered_ipd_sat
centered_colnames <- grep("_CENTERED$", colnames(centered_ipd_sat), value = TRUE)
weighted_data <- estimate_weights(data = centered_ipd_sat, centered_colnames = centered_colnames)
weighted_data2 <- estimate_weights(
  data = centered_ipd_sat, centered_colnames = centered_colnames,
  n_boot_iteration = 100
)

# get dummy binary IPD
pseudo_adrs <- get_pseudo_ipd_binary(
  binary_agd = data.frame(
    ARM = rep("B", 2),
    RESPONSE = c("YES", "NO"),
    COUNT = c(280, 120)
  ),
  format = "stacked"
)

# unanchored binary MAIC, with CI based on sandwich estimator
maic_unanchored(
  weights_object = weighted_data,
  ipd = adrs_sat,
  pseudo_ipd = pseudo_adrs,
  trt_ipd = "A",
  trt_agd = "B",
  trt_var_ipd = "ARM",
  trt_var_agd = "ARM",
  endpoint_type = "binary",
  endpoint_name = "Binary Endpoint",
  eff_measure = "RR",
  # binary specific args
  binary_robust_cov_type = "HC3"
)

# unanchored binary MAIC, with bootstrapped CI
maic_unanchored(
  weights_object = weighted_data2,
  ipd = adrs_sat,
  pseudo_ipd = pseudo_adrs,
  trt_ipd = "A",
  trt_agd = "B",
  trt_var_ipd = "ARM",
  trt_var_agd = "ARM",
  endpoint_type = "binary",
  endpoint_name = "Binary Endpoint",
  eff_measure = "RR",
  # binary specific args
  binary_robust_cov_type = "HC3"
)

#---------------------------------
