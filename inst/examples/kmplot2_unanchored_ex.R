# unanchored example using kmplot2
data(weighted_sat)
data(adtte_sat)
data(pseudo_ipd_sat)

kmplot2(
  weights_object = weighted_sat,
  tte_ipd = adtte_sat,
  tte_pseudo_ipd = pseudo_ipd_sat,
  trt_ipd = "A",
  trt_agd = "B",
  trt_common = NULL,
  trt_var_ipd = "ARM",
  trt_var_agd = "ARM",
  endpoint_name = "Overall Survival",
  km_conf_type = "log-log",
  time_scale = "month",
  break_x_by = 2,
  xlim = c(0, 20)
)
