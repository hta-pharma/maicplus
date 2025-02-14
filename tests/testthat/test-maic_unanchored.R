test_that("test binary case", {
  data(centered_ipd_sat)
  data(agd)
  agd <- process_agd(agd)

  ipd_centered <- center_ipd(ipd = centered_ipd_sat, agd = agd)

  # estimate weights
  centered_colnames <- c("AGE", "AGE_SQUARED", "SEX_MALE", "ECOG0", "SMOKE", "N_PR_THER_MEDIAN")
  centered_colnames <- paste0(centered_colnames, "_CENTERED")

  weighted_data <- estimate_weights(data = ipd_centered, centered_colnames = centered_colnames)
  weighted_data2 <- estimate_weights(
    data = ipd_centered, centered_colnames = centered_colnames, n_boot_iteration = 20,
    set_seed_boot = 1234
  )
  # get dummy binary IPD
  data(adrs_sat)

  pseudo_adrs <- get_pseudo_ipd_binary(
    binary_agd = data.frame(
      ARM = rep("B", 2),
      RESPONSE = c("YES", "NO"),
      COUNT = c(280, 120)
    ),
    format = "stacked"
  )

  # unanchored binary MAIC, with CI based on sandwich estimator
  testout_RR <-
    maic_unanchored(
      weights_object = weighted_data,
      ipd = adrs_sat,
      pseudo_ipd = pseudo_adrs,
      trt_ipd = "A",
      trt_agd = "B",
      trt_var_ipd = "ARM",
      trt_var_agd = "ARM",
      endpoint_type = "binary",
      endpoint_name = "Binary Endpoint",
      eff_measure = "RR",
      # binary specific args
      binary_robust_cov_type = "HC3"
    )

  testout_RD <-
    maic_unanchored(
      weights_object = weighted_data,
      ipd = adrs_sat,
      pseudo_ipd = pseudo_adrs,
      trt_ipd = "A",
      trt_agd = "B",
      trt_var_ipd = "ARM",
      trt_var_agd = "ARM",
      endpoint_type = "binary",
      endpoint_name = "Binary Endpoint",
      eff_measure = "RD",
      # binary specific args
      binary_robust_cov_type = "HC3"
    )


  testout_OR <-
    maic_unanchored(
      weights_object = weighted_data,
      ipd = adrs_sat,
      pseudo_ipd = pseudo_adrs,
      trt_ipd = "A",
      trt_agd = "B",
      trt_var_ipd = "ARM",
      trt_var_agd = "ARM",
      endpoint_type = "binary",
      endpoint_name = "Binary Endpoint",
      eff_measure = "OR",
      # binary specific args
      binary_robust_cov_type = "HC3"
    )


  # unanchored binary MAIC, with bootstrapped CI
  testout_boot_RR <-
    maic_unanchored(
      weights_object = weighted_data2,
      ipd = adrs_sat,
      pseudo_ipd = pseudo_adrs,
      trt_ipd = "A",
      trt_agd = "B",
      trt_var_ipd = "ARM",
      trt_var_agd = "ARM",
      endpoint_type = "binary",
      endpoint_name = "Binary Endpoint",
      eff_measure = "RR",
      # binary specific args
      binary_robust_cov_type = "HC3"
    )


  # Compare robust outputs
  expect_snapshot(testout_RR$descriptive$summary)
  expect_snapshot(testout_RR$inferential$summary)
  expect_snapshot(testout_RR$inferential$fit)

  expect_snapshot(testout_RD$descriptive$summary)
  expect_snapshot(testout_RD$inferential$summary)
  expect_snapshot(testout_RD$inferential$fit)

  expect_snapshot(testout_OR$descriptive$summary)
  expect_snapshot(testout_OR$inferential$summary)
  expect_snapshot(testout_OR$inferential$fit)

  # Compare bootstrap outputs
  expect_snapshot(print(testout_boot_RR$descriptive$summary, digits = 5))
  expect_snapshot(print(testout_boot_RR$inferential$summary, digits = 5))
  expect_snapshot(print(testout_boot_RR$inferential$fit, digits = 5))
})



test_that("check if match properly", {
  data(centered_ipd_twt)
  data(adrs_twt)

  centered_colnames <- c("AGE", "AGE_SQUARED", "SEX_MALE", "ECOG0", "SMOKE", "N_PR_THER_MEDIAN")
  centered_colnames <- paste0(centered_colnames, "_CENTERED")

  weighted_data <- estimate_weights(
    data = centered_ipd_twt,
    centered_colnames = centered_colnames
  )

  # get dummy binary IPD
  pseudo_adrs <- get_pseudo_ipd_binary(
    binary_agd = data.frame(
      ARM = c("B", "C", "B", "C"),
      RESPONSE = c("YES", "YES", "NO", "NO"),
      COUNT = c(280, 120, 200, 200)
    ),
    format = "stacked"
  )

  testout2 <- maic_anchored(
    weights_object = weighted_data,
    ipd = adrs_twt,
    pseudo_ipd = pseudo_adrs,
    trt_ipd = "A",
    trt_agd = "B",
    trt_common = "C",
    normalize_weight = FALSE,
    endpoint_type = "binary",
    endpoint_name = "Binary Endpoint",
    eff_measure = "OR",
    # binary specific args
    binary_robust_cov_type = "HC3"
  )

  weights_before_wrapper <- weighted_data$data$weights
  adrs_twt_dummy <- adrs_twt
  adrs_twt_dummy$weights <- weighted_data$data$weights[match(adrs_twt_dummy$USUBJID, weighted_data$data$USUBJID)]

  weights_after_wrapper <- testout2$inferential$fit$model_after$data$weights
  expect_equal(adrs_twt_dummy$weights, weights_after_wrapper)
})

test_that("test time to event case", {
  data(centered_ipd_sat)
  data(agd)
  agd <- process_agd(agd)

  ipd_centered <- center_ipd(ipd = centered_ipd_sat, agd = agd)

  # estimate weights
  centered_colnames <- c("AGE", "AGE_SQUARED", "SEX_MALE", "ECOG0", "SMOKE", "N_PR_THER_MEDIAN")
  centered_colnames <- paste0(centered_colnames, "_CENTERED")

  weighted_data <- estimate_weights(
    data = ipd_centered,
    centered_colnames = centered_colnames,
    start_val = 0,
    method = "BFGS"
  )

  weighted_data_boot <- estimate_weights(
    data = ipd_centered,
    centered_colnames = centered_colnames,
    start_val = 0,
    method = "BFGS",
    n_boot_iteration = 500,
    set_seed_boot = 1234
  )

  # inferential result
  testout <- maic_unanchored(
    weights_object = weighted_data,
    ipd = adtte_sat,
    pseudo_ipd = pseudo_ipd_sat,
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

  testout2 <- maic_unanchored(
    weights_object = weighted_data_boot,
    ipd = adtte_sat,
    pseudo_ipd = pseudo_ipd_sat,
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


  # Compare robust outputs
  expect_snapshot(testout$descriptive$summary)
  expect_snapshot(testout$inferential$summary)
  expect_snapshot(testout$inferential$fit)

  # Compare bootstrap outputs
  expect_snapshot(print(testout2$descriptive$summary, digits = 5))
  expect_snapshot(print(testout2$inferential$summary, digits = 5))
  expect_snapshot(print(testout2$inferential$fit, digits = 5))
})
