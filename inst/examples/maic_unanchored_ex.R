#
# unanchored example using maic_unanchored for time-to-event data
#
data(centered_ipd_sat)
data(adtte_sat)
data(pseudo_ipd_sat)

#### derive weights
weighted_data <- estimate_weights(
  data = centered_ipd_sat,
  centered_colnames = grep("_CENTERED$", names(centered_ipd_sat)),
  start_val = 0,
  method = "BFGS"
)

weighted_data2 <- estimate_weights(
  data = centered_ipd_sat,
  centered_colnames = grep("_CENTERED$", names(centered_ipd_sat)),
  start_val = 0,
  method = "BFGS",
  n_boot_iteration = 100,
  set_seed_boot = 1234
)

# inferential result
result <- maic_unanchored(
  weights_object = weighted_data,
  ipd = adtte_sat,
  pseudo_ipd = pseudo_ipd_sat,
  trt_var_ipd = "ARM",
  trt_var_agd = "ARM",
  trt_ipd = "A",
  trt_agd = "B",
  endpoint_name = "Overall Survival",
  endpoint_type = "tte",
  eff_measure = "HR",
  time_scale = "month",
  km_conf_type = "log-log"
)
result$descriptive$summary
result$inferential$summary

result_boot <- maic_unanchored(
  weights_object = weighted_data2,
  ipd = adtte_sat,
  pseudo_ipd = pseudo_ipd_sat,
  trt_var_ipd = "ARM",
  trt_var_agd = "ARM",
  trt_ipd = "A",
  trt_agd = "B",
  endpoint_name = "Overall Survival",
  endpoint_type = "tte",
  eff_measure = "HR",
  time_scale = "month",
  km_conf_type = "log-log"
)
result$descriptive$summary
result$inferential$summary
