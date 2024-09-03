# Anchored example using maic_anchored for binary outcome
data(weighted_twt)
data(adrs_twt)

# Reported summary data
pseudo_adrs <- get_pseudo_ipd_binary(
  binary_agd = data.frame(
    ARM = c("B", "C", "B", "C"),
    RESPONSE = c("YES", "YES", "NO", "NO"),
    COUNT = c(280, 120, 200, 200)
  ),
  format = "stacked"
)

# inferential result
result_binary <- maic_anchored(
  weights_object = weighted_twt,
  ipd = adrs_twt,
  pseudo_ipd = pseudo_adrs,
  trt_var_ipd = "ARM",
  trt_var_agd = "ARM",
  trt_ipd = "A",
  trt_agd = "B",
  trt_common = "C",
  endpoint_name = "Binary Event",
  endpoint_type = "binary",
  eff_measure = "OR"
)

result_binary$descriptive$summary
result_binary$inferential$summary
