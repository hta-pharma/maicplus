# anchored example using kmplot
data(weighted_twt)
data(adtte_twt)
data(pseudo_ipd_twt)

# plot by trial
kmplot(
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
  time_grid = seq(0, 20, by = 2),
  use_colors = NULL,
  use_line_types = NULL,
  use_pch_cex = 0.65,
  use_pch_alpha = 100
)

# plot by arm
kmplot(
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
  time_grid = seq(0, 20, by = 2),
  use_colors = NULL,
  use_line_types = NULL,
  use_pch_cex = 0.65,
  use_pch_alpha = 100
)

# plot all
kmplot(
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
  time_grid = seq(0, 20, by = 2),
  use_colors = NULL,
  use_line_types = NULL,
  use_pch_cex = 0.65,
  use_pch_alpha = 100
)
