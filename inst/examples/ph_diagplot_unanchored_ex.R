# unanchored example using ph_diagplot
data(weighted_sat)
data(adtte_sat)
data(pseudo_ipd_sat)

ph_diagplot(
  weights_object = weighted_sat,
  tte_ipd = adtte_sat,
  tte_pseudo_ipd = pseudo_ipd_sat,
  trt_var_ipd = "ARM",
  trt_var_agd = "ARM",
  trt_ipd = "A",
  trt_agd = "B",
  trt_common = NULL,
  endpoint_name = "Overall Survival",
  time_scale = "week",
  zph_transform = "log",
  zph_log_hazard = TRUE
)
