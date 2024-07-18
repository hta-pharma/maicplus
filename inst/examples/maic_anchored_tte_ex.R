# Anchored example using maic_anchored for time-to-event data

weighted_twt # previously computed result from estimate_weights()
adtte_twt # two arm time-to-event data

pseudo_ipd_twt # digitised time-to-event data

result <- maic_anchored(
  weights_object = weighted_twt,
  ipd = adtte_twt,
  trt_var_ipd = "ARM",
  pseudo_ipd = pseudo_ipd_twt,
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
result$inferential$report_median_surv
result$inferential$report_overall_robustCI
result$inferential$report_overall_bootCI
