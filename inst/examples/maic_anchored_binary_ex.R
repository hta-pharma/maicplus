# anchored example using maic_anchored for binary outcome

# Uses example data
centered_ipd_twt # patient characteristics data
adrs_twt # outcome data

# get dummy binary IPD
pseudo_adrs <- get_pseudo_ipd_binary(
  binary_agd = data.frame(
    ARM = c("B", "C", "B", "C"),
    RESPONSE = c("YES", "YES", "NO", "NO"),
    COUNT = c(280, 120, 200, 200)
  ),
  format = "stacked"
)

#### derive weights
match_res <- estimate_weights(
  data = centered_ipd_twt,
  centered_colnames = grep("_CENTERED$", names(centered_ipd_twt)),
  start_val = 0,
  method = "BFGS",
  n_boot_iteration = 10,
  set_seed_boot = 1234
)

# inferential result
result <- maic_anchored(
  weights_object = match_res,
  ipd = adrs_twt,
  trt_var_ipd = "ARM",
  pseudo_ipd = pseudo_adrs,
  trt_var_agd = "ARM",
  trt_ipd = "A",
  trt_agd = "B",
  trt_common = "C",
  endpoint_name = "Binary Event",
  endpoint_type = "binary",
  eff_measure = "OR"
)

result$inferential$report_overall_robustCI
result$inferential$report_overall_bootCI
