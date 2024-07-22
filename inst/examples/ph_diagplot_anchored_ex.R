# anchored example using ph_diagplot
data(weighted_twt)
data(adtte_twt)
data(pseudo_ipd_twt)

ph_diagplot(
  weights_object = weighted_twt,
  tte_ipd = adtte_twt,
  tte_pseudo_ipd = pseudo_ipd_twt,
  trt_var_ipd = "ARM",
  trt_var_agd = "ARM",
  trt_ipd = "A",
  trt_agd = "B",
  trt_common = "C",
  endpoint_name = "Overall Survival",
  time_scale = "week",
  zph_transform = "log",
  zph_log_hazard = TRUE
)
