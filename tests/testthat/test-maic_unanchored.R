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
  testout <-
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

  # unanchored binary MAIC, with bootstrapped CI
  testout2 <-
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

  if (FALSE) {
    # Manual snapshot of results
    expectout <- testout
    expectout2 <- testout2
    save(list = c("expectout", "expectout2"), file = test_path("data", "test_binary_unanchored_expected.RData"))
  }

  load(test_path("data", "test_binary_unanchored_expected.RData"))

  # Compare robust outputs
  expect_equal(testout$descriptive, expectout$descriptive)
  expect_equal(testout$inferential, expectout$inferential)
  expect_equal(testout$inferential$model_before, expectout$inferential$model_before)
  expect_equal(testout$inferential$model_after, expectout$inferential$model_after)

  # Compare bootstrap outputs
  expect_equal(testout2$descriptive, expectout2$descriptive)
  expect_equal(testout2$inferential$boot_est["t"], expectout2$inferential$boot_est["t"])
  expect_equal(testout2$inferential$boot_est["seed"], expectout2$inferential$boot_est["seed"])
  expect_equal(testout2$inferential$report_overall_bootCI, expectout2$inferential$report_overall_bootCI)
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
  
  if (FALSE) {
    # Manual snapshot of results
    expectout <- testout
    expectout2 <- testout2
    save(list = c("expectout", "expectout2"), file = test_path("data", "test_tte_unanchored_expected.RData"))
  }
  
  load(test_path("data", "test_tte_unanchored_expected.RData"))
  # Compare robust outputs
  expect_equal(testout$descriptive, expectout$descriptive)
  expect_equal(testout$inferential, expectout$inferential)
  expect_equal(testout$inferential$model_before, expectout$inferential$model_before)
  expect_equal(testout$inferential$model_after, expectout$inferential$model_after)
  
  # Compare bootstrap outputs
  expect_equal(testout2$descriptive, expectout2$descriptive)
  expect_equal(testout2$inferential$boot_est["t"], expectout2$inferential$boot_est["t"])
  expect_equal(testout2$inferential$boot_est["seed"], expectout2$inferential$boot_est["seed"])
  expect_equal(testout2$inferential$report_overall_bootCI, expectout2$inferential$report_overall_bootCI)
})
