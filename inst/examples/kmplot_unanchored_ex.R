# unanchored example using kmplot
data(weighted_sat)
data(adtte_sat)
data(pseudo_ipd_sat)

kmplot(
  weights_object = weighted_sat,
  tte_ipd = adtte_sat,
  tte_pseudo_ipd = pseudo_ipd_sat,
  trt_var_ipd = "ARM",
  trt_var_agd = "ARM",
  endpoint_name = "Overall Survival",
  trt_ipd = "A",
  trt_agd = "B",
  trt_common = NULL,
  km_conf_type = "log-log",
  time_scale = "month",
  time_grid = seq(0, 20, by = 2),
  use_colors = NULL,
  use_line_types = NULL,
  use_pch_cex = 0.65,
  use_pch_alpha = 100
)
