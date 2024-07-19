# Anchored example using maic_anchored for binary outcome

# weighted_twt    # weighted patient level data from estimate_weights()
# adrs_twt        # patient level outcome data
# pseudo_adrs     # pseudo patient level outcome data from publication

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

result_binary$inferential$report_overall_robustCI
result_binary$inferential$report_overall_bootCI


# Anchored example using maic_anchored for time-to-event data

# weighted_twt    # weighted patient level data from estimate_weights()
# adtte_twt       # two arm patient level time-to-event data
# pseudo_ipd_twt  # digitised time-to-event data

result_tte <- maic_anchored(
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
result_tte$inferential$report_median_surv
result_tte$inferential$report_overall_robustCI
result_tte$inferential$report_overall_bootCI
