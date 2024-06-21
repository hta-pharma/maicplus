# anchored example using kmplot

### IPD
# Read in relevant ADaM data and rename variables of interest
adsl_twt
adtte_twt

### AgD
# Baseline aggregate data for the comparator population
agd

# for time-to-event endpoints, pseudo IPD from digitalized KM
pseudo_ipd_twt

#### prepare data
target_pop <- process_agd(agd)
adsl <- dummize_ipd(adsl_twt, dummize_cols = c("SEX"), dummize_ref_level = c("Female"))
use_adsl <- center_ipd(ipd = adsl, agd = target_pop)

#### derive weights
match_res <- estimate_weights(
  data = use_adsl,
  centered_colnames = grep("_CENTERED$", names(use_adsl)),
  start_val = 0,
  method = "BFGS"
)

# plot by trial
kmplot(
  weights_object = match_res,
  tte_ipd = adtte_twt,
  trt_var_ipd = "ARM",
  tte_pseudo_ipd = pseudo_ipd_twt,
  trt_var_agd = "ARM",
  endpoint_name = "Overall Survival",
  trt_ipd = "A",
  trt_agd = "B",
  trt_common = "C",
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
  weights_object = match_res,
  tte_ipd = adtte_twt,
  trt_var_ipd = "ARM",
  tte_pseudo_ipd = pseudo_ipd_twt,
  trt_var_agd = "ARM",
  endpoint_name = "Overall Survival",
  trt_ipd = "A",
  trt_agd = "B",
  trt_common = "C",
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
  weights_object = match_res,
  tte_ipd = adtte_twt,
  trt_var_ipd = "ARM",
  tte_pseudo_ipd = pseudo_ipd_twt,
  trt_var_agd = "ARM",
  endpoint_name = "Overall Survival",
  trt_ipd = "A",
  trt_agd = "B",
  trt_common = "C",
  km_conf_type = "log-log",
  km_layout = "all",
  time_scale = "month",
  time_grid = seq(0, 20, by = 2),
  use_colors = NULL,
  use_line_types = NULL,
  use_pch_cex = 0.65,
  use_pch_alpha = 100
)
