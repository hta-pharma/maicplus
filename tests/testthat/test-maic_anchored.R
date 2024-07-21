test_that("maic_anchored works for TTE using robust SE", {
  
  data(adtte_twt)
  data(pseudo_ipd_twt)
  data(centered_ipd_twt)
  
  #### derive weights
  cols <- grep("_CENTERED$", names(centered_ipd_twt), value = TRUE)
  weighted_data <- estimate_weights(
    data = centered_ipd_twt,
    centered_colnames = cols,
    start_val = 0,
    method = "BFGS"
  )
  
  # inferential result
  result <- maic_anchored(
    weights_object = weighted_data,
    ipd = adtte_twt,
    pseudo_ipd = pseudo_ipd_twt,
    trt_var_ipd = "ARM",
    trt_var_agd = "ARM",
    trt_ipd = "A",
    trt_agd = "B",
    trt_common = "C",
    endpoint_name = "Overall Survival",
    endpoint_type = "tte",
    eff_measure = "HR",
    time_scale = "month",
    km_conf_type = "log-log"
  )
  
  testthat::expect_true(is.list(result["descriptive"]))
  testthat::expect_true(is.list(result["inferential"]))
  expect_equal(
    result$inferential$report_median_surv$rmean,
    c(
      2.5647965487863, 8.7096897111058, 2.6794733007394, 10.5846090795552,
      2.4552717139068, 4.3035505695334
    )
  )
  expect_equal(
    result$inferential$report_median_surv$`se(rmean)`,
    c(
      0.113669935856188, 0.355147660158620, 0.206708266956460, 0.573979370096874,
      0.098488879305725, 0.336726020204787
    )
  )
  expect_equal(
    result$inferential$report_overall_robustCI$`median[95% CI]`,
    c("7.6[6.3;10.3]", "1.8[1.6; 2.0]", "12.2[10.2; NA]", " 1.8[ 1.5;2.3]", "2.7[2.3;3.3]", "1.9[1.7;2.0]", "--")
  )
  expect_equal(
    result$inferential$report_overall_robustCI$`HR[95% CI]`,
    c("0.22[0.19;0.26]", "", "0.16[0.11;0.24]", "", "0.57[0.48;0.68]", "", "0.29 [0.19; 0.43]")
  )
})
