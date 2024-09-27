# Anchored example using maic_anchored for time-to-event data
data(weighted_twt)
data(adtte_twt)
data(pseudo_ipd_twt)

result_tte <- maic_anchored(
  weights_object = weighted_twt,
  ipd = adtte_twt,
  pseudo_ipd = pseudo_ipd_twt,
  trt_var_ipd = "ARM",
  trt_var_agd = "ARM",
  trt_ipd = "A",
  trt_agd = "B",
  trt_common = "C",
  endpoint_name = "Overall Survival",
  endpoint_type = "tte",
  eff_measure = "HR",
  time_scale = "month",
  km_conf_type = "log-log",
)
result_tte$descriptive$summary
result_tte$inferential$summary
