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

test_that("maic_anchored works for TTE using bootstrap SE", {
  data(adtte_twt)
  data(pseudo_ipd_twt)
  data(centered_ipd_twt)

  #### derive weights
  cols <- grep("_CENTERED$", names(centered_ipd_twt), value = TRUE)
  weighted_data_boot <- estimate_weights(
    data = centered_ipd_twt,
    centered_colnames = cols,
    start_val = 0,
    method = "BFGS",
    n_boot_iteration = 5,
    set_seed_boot = 1234
  )

  result <- maic_anchored(
    weights_object = weighted_data_boot,
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

  expect_equal(
    result$inferential$report_overall_bootCI$`median[95% CI]`,
    c("7.6[6.3;10.3]", "1.8[1.6; 2.0]", "12.2[10.2; NA]", " 1.8[ 1.5;2.3]", "2.7[2.3;3.3]", "1.9[1.7;2.0]", "--")
  )
  expect_equal(
    result$inferential$report_overall_bootCI$`HR[95% CI]`,
    c("0.22[0.19;0.26]", "", "0.16[0.11;0.24]", "", "0.57[0.48;0.68]", "", "0.29 [0.21; 0.55]")
  )

  t_matrix_expected <- matrix(
    c(
      -1.5224191255192, 0.060103109952175, 0.24515935624033, -2.0813843502016, 0.22881474544327, 0.052356187732271,
      -1.0124623878469, 0.040408249516049, 0.20101803281310, -1.5714276125293, 0.18072445129573, 0.032661327296144,
      -1.6841548337608, 0.056582360570672, 0.23787047015271, -2.2431200584432, 0.22098741672495, 0.048835438350767,
      -1.4861313293074, 0.045824723120669, 0.21406709957550, -2.0450965539898, 0.19513533995862, 0.038077800900765,
      -1.3989135763837, 0.035646484583158, 0.18880276635462, -1.9578788010661, 0.16703162084843, 0.027899562363253
    ),
    byrow = TRUE,
    ncol = 6
  )
  expect_equal(result$inferential$boot_est$t, t_matrix_expected)
})


test_that("maic_anchored for binary case gives the expected result", {
  data(weighted_twt)
  data(adrs_twt)

  # Reported summary data
  pseudo_adrs <- get_pseudo_ipd_binary(
    binary_agd = data.frame(
      ARM = c("B", "C", "B", "C"),
      RESPONSE = c("YES", "YES", "NO", "NO"),
      COUNT = c(280, 120, 200, 200)
    ),
    format = "stacked"
  )

  # inferential result
  expect_warning(
    result <- maic_anchored(
      weights_object = weighted_twt,
      ipd = adrs_twt,
      pseudo_ipd = pseudo_adrs,
      trt_var_ipd = "ARM",
      trt_var_agd = "ARM",
      trt_ipd = "A",
      trt_agd = "B",
      trt_common = "C",
      endpoint_name = "Binary Event",
      endpoint_type = "binary",
      eff_measure = "OR"
    ),
    "non-integer"
  )

  expect_equal(
    result$inferential$report_overall_robustCI$`OR[95% CI]`,
    c("1.70[1.28;2.26]", "", "1.14[0.67;1.95]", "", "2.33[1.75;3.12]", "", "0.49 [0.27; 0.90]")
  )
  expect_equal(
    result$inferential$report_overall_robustCI$`n.events(%)`,
    c("390(78.0)", "338(67.6)", "128.8(25.8)", "124.2(24.8)", "280(58.3)", "120(37.5)", "--")
  )
  expect_equal(
    result$inferential$report_overall_robustCI$N,
    c("500", "500", "500", "500", "480", "320", "--")
  )
  expect_equal(
    result$inferential$report_overall_bootCI$`OR[95% CI]`,
    c("1.70[1.28;2.26]", "", "1.14[0.33;0.98]", "", "2.33[1.75;3.12]", "", "0.49 [0.14; 0.42]")
  )
})
