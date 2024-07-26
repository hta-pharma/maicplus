# anchored example using kmplot2
data(weighted_twt)
data(adtte_twt)
data(pseudo_ipd_twt)

# plot by trial
kmplot2(
  weights_object = weighted_twt,
  tte_ipd = adtte_twt,
  tte_pseudo_ipd = pseudo_ipd_twt,
  trt_ipd = "A",
  trt_agd = "B",
  trt_common = "C",
  trt_var_ipd = "ARM",
  trt_var_agd = "ARM",
  endpoint_name = "Overall Survival",
  km_conf_type = "log-log",
  km_layout = "by_trial",
  time_scale = "month",
  break_x_by = 2
)

# plot by arm
kmplot2(
  weights_object = weighted_twt,
  tte_ipd = adtte_twt,
  tte_pseudo_ipd = pseudo_ipd_twt,
  trt_ipd = "A",
  trt_agd = "B",
  trt_common = "C",
  trt_var_ipd = "ARM",
  trt_var_agd = "ARM",
  endpoint_name = "Overall Survival",
  km_conf_type = "log-log",
  km_layout = "by_arm",
  time_scale = "month",
  break_x_by = 2
)

# plot all
kmplot2(
  weights_object = weighted_twt,
  tte_ipd = adtte_twt,
  tte_pseudo_ipd = pseudo_ipd_twt,
  trt_ipd = "A",
  trt_agd = "B",
  trt_common = "C",
  trt_var_ipd = "ARM",
  trt_var_agd = "ARM",
  endpoint_name = "Overall Survival",
  km_conf_type = "log-log",
  km_layout = "all",
  time_scale = "month",
  break_x_by = 2,
  xlim = c(0, 20),
  show_risk_set = FALSE
)
