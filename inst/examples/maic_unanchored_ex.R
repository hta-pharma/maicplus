# anchored example using maic_anchored for tte



# Read in relevant ADaM data and rename variables of interest
adsl <- read.csv(system.file("extdata", "adsl.csv",
  package = "maicplus",
  mustWork = TRUE
))
adtte <- read.csv(system.file("extdata", "adtte.csv",
  package = "maicplus",
  mustWork = TRUE
))
adtte$TIME <- adtte$AVAL
adtte$EVENT <- adtte$EVNT
adtte <- adtte[adtte$ARM == "A", , drop = FALSE]
adsl <- adsl[adsl$USUBJID %in% adtte$USUBJID, , drop = FALSE]

### AgD
# Baseline aggregate data for the comparator population
target_pop <- read.csv(system.file("extdata", "aggregate_data_example_1.csv",
  package = "maicplus", mustWork = TRUE
))
# for time-to-event endpoints, pseudo IPD from digitalized KM
pseudo_ipd <- read.csv(system.file("extdata", "psuedo_IPD.csv",
  package = "maicplus",
  mustWork = TRUE
))
pseudo_ipd$ARM <- "B"

#### prepare data
target_pop <- process_agd(target_pop)
adsl <- dummize_ipd(adsl, dummize_cols = c("SEX"), dummize_ref_level = c("Female"))
use_adsl <- center_ipd(ipd = adsl, agd = target_pop)

#### derive weights
match_res <- estimate_weights(
  data = use_adsl,
  centered_colnames = grep("_CENTERED$", names(use_adsl)),
  start_val = 0,
  method = "BFGS"
)

match_res_boot <- estimate_weights(
  data = use_adsl,
  centered_colnames = grep("_CENTERED$", names(use_adsl)),
  start_val = 0,
  method = "BFGS",
  n_boot_iteration = 500,
  set_seed_boot = 1234
)

# inferential result
result <- maic_unanchored(
  weights_object = match_res,
  ipd = adtte,
  pseudo_ipd = pseudo_ipd,
  trt_var_ipd = "ARM",
  trt_var_agd = "ARM",
  trt_ipd = "A",
  trt_agd = "B",
  endpoint_name = "Overall Survival",
  endpoint_type = "tte",
  eff_measure = "HR",
  time_scale = "month",
  km_conf_type = "log-log"
)
result$inferential$report_median_surv
result$inferential$report_overall_robustCI
result$inferential$report_overall_bootCI

result_boot <- maic_unanchored(
  weights_object = match_res_boot,
  ipd = adtte,
  pseudo_ipd = pseudo_ipd,
  trt_var_ipd = "ARM",
  trt_var_agd = "ARM",
  trt_ipd = "A",
  trt_agd = "B",
  endpoint_name = "Overall Survival",
  endpoint_type = "tte",
  eff_measure = "HR",
  time_scale = "month",
  km_conf_type = "log-log"
)
result_boot$inferential$report_median_surv
result_boot$inferential$report_overall_robustCI
result_boot$inferential$report_overall_bootCI
quantile(result_boot$inferential$boot_est, p = 0.025)
quantile(result_boot$inferential$boot_est, p = 0.975)

ph_diagplot(
  weights_object = match_res,
  tte_ipd = adtte,
  tte_pseudo_ipd = pseudo_ipd,
  trt_var_ipd = "ARM",
  trt_var_agd = "ARM",
  trt_ipd = "A",
  trt_agd = "B",
  endpoint_name = "Overall Survival",
  time_scale = "week",
  zph_transform = "log",
  zph_log_hazard = TRUE
)

kmplot(
  weights_object = match_res,
  tte_ipd = adtte,
  tte_pseudo_ipd = pseudo_ipd,
  trt_var_ipd = "ARM",
  trt_var_agd = "ARM",
  endpoint_name = "Overall Survival",
  trt_ipd = "A",
  trt_agd = "B",
  km_conf_type = "log-log",
  km_layout = "by_trial",
  time_scale = "month",
  time_grid = seq(0, 20, by = 2),
  use_colors = NULL,
  use_line_types = NULL,
  use_pch_cex = 0.65,
  use_pch_alpha = 100
)
