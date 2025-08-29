

#setwd("~/GitHub/maicplus")
source("design/explore_normalization/normalize_weights.R")

library(maicplus)
library(dplyr)

data(centered_ipd_sat)
data(adtte_sat)
data(pseudo_ipd_sat)

centered_colnames <- c("AGE", "AGE_SQUARED", "SEX_MALE", "ECOG0", "SMOKE", "N_PR_THER_MEDIAN")
centered_colnames <- paste0(centered_colnames, "_CENTERED")

weighted_data <- estimate_weights(
  data = centered_ipd_sat,
  centered_colnames = centered_colnames
)
weighted_data_SW1 <- weighted_data_SWN <- weighted_data_SWESS <- weighted_data
weighted_data_SW1$data$weights <- normalize_weights(weighted_data$data$weights, method = "SW1")
weighted_data_SWN$data$weights <- normalize_weights(weighted_data$data$weights, method = "SWN")
weighted_data_SWESS$data$weights <- normalize_weights(weighted_data$data$weights, method = "SWESS")


### OW
result_OW <- maic_unanchored(
  weights_object = weighted_data,
  ipd = adtte_sat,
  pseudo_ipd = pseudo_ipd_sat,
  trt_ipd = "A",
  trt_agd = "B",
  normalize_weight = FALSE,
  endpoint_name = "Overall Survival",
  endpoint_type = "tte",
  eff_measure = "HR",
  time_scale = "month",
  km_conf_type = "log-log"
)

### SW1
result_SW1 <- maic_unanchored(
  weights_object = weighted_data_SW1,
  ipd = adtte_sat,
  pseudo_ipd = pseudo_ipd_sat,
  trt_ipd = "A",
  trt_agd = "B",
  normalize_weight = FALSE,
  endpoint_name = "Overall Survival",
  endpoint_type = "tte",
  eff_measure = "HR",
  time_scale = "month",
  km_conf_type = "log-log"
)

### SWN
result_SWN <- maic_unanchored(
  weights_object = weighted_data_SWN,
  ipd = adtte_sat,
  pseudo_ipd = pseudo_ipd_sat,
  trt_ipd = "A",
  trt_agd = "B",
  normalize_weight = FALSE,
  endpoint_name = "Overall Survival",
  endpoint_type = "tte",
  eff_measure = "HR",
  time_scale = "month",
  km_conf_type = "log-log"
)

### SWESS
result_SWESS <- maic_unanchored(
  weights_object = weighted_data_SWESS,
  ipd = adtte_sat,
  pseudo_ipd = pseudo_ipd_sat,
  trt_ipd = "A",
  trt_agd = "B",
  normalize_weight = FALSE,
  endpoint_name = "Overall Survival",
  endpoint_type = "tte",
  eff_measure = "HR",
  time_scale = "month",
  km_conf_type = "log-log"
)

# results are not the same
result_OW$inferential$summary
result_SW1$inferential$summary
result_SWN$inferential$summary
result_SWESS$inferential$summary

# baseline hazards are different
library(survival)
basehaz(result_OW$inferential$fit$model_after)[1:10,]
basehaz(result_SW1$inferential$fit$model_after)[1:10,]

basehaz(result_OW$inferential$fit$model_before)[1:10,]
basehaz(result_SW1$inferential$fit$model_before)[1:10,]
