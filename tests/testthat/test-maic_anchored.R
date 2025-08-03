test_that("maic_anchored works for TTE", {
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

  weighted_data_boot <- estimate_weights(
    data = centered_ipd_twt,
    centered_colnames = cols,
    start_val = 0,
    method = "BFGS",
    n_boot_iteration = 5,
    set_seed_boot = 1234
  )

  # inferential result
  testout <- maic_anchored(
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

  testout2 <- maic_anchored(
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

  # Compare robust outputs
  expect_snapshot(testout$descriptive$summary)
  expect_snapshot(print(testout$inferential$summary, digits = 5))
  expect_snapshot(testout$inferential$fit)

  # Compare bootstrap outputs
  expect_snapshot(print(testout2$descriptive$summary, digits = 5))
  expect_snapshot(print(testout2$inferential$summary, digits = 5))
  expect_snapshot(print(testout2$inferential$fit, digits = 5))
})


test_that("maic_anchored for binary case gives the expected result", {
  data(centered_ipd_twt)
  data(agd)
  agd <- process_agd(agd)

  ipd_centered <- center_ipd(ipd = centered_ipd_twt, agd = agd)

  # estimate weights
  centered_colnames <- c("AGE", "AGE_SQUARED", "SEX_MALE", "ECOG0", "SMOKE", "N_PR_THER_MEDIAN")
  centered_colnames <- paste0(centered_colnames, "_CENTERED")

  weighted_data <- estimate_weights(data = ipd_centered, centered_colnames = centered_colnames)
  weighted_data2 <- estimate_weights(
    data = ipd_centered, centered_colnames = centered_colnames,
    n_boot_iteration = 20, set_seed_boot = 1234
  )
  # get dummy binary IPD
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
  testout_OR <- maic_anchored(
    weights_object = weighted_data,
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
  )

  testout_RR <- maic_anchored(
    weights_object = weighted_data,
    ipd = adrs_twt,
    pseudo_ipd = pseudo_adrs,
    trt_var_ipd = "ARM",
    trt_var_agd = "ARM",
    trt_ipd = "A",
    trt_agd = "B",
    trt_common = "C",
    endpoint_name = "Binary Event",
    endpoint_type = "binary",
    eff_measure = "RR"
  )

  testout_RD <- maic_anchored(
    weights_object = weighted_data,
    ipd = adrs_twt,
    pseudo_ipd = pseudo_adrs,
    trt_var_ipd = "ARM",
    trt_var_agd = "ARM",
    trt_ipd = "A",
    trt_agd = "B",
    trt_common = "C",
    endpoint_name = "Binary Event",
    endpoint_type = "binary",
    eff_measure = "RD"
  )

  # bootstrap
  testout_boot_OR <- maic_anchored(
    weights_object = weighted_data2,
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
  )

  # Compare robust outputs
  expect_snapshot(testout_OR$descriptive$summary)
  expect_snapshot(testout_OR$inferential$summary)
  expect_snapshot(testout_OR$inferential$fit)

  expect_snapshot(testout_RR$descriptive$summary)
  expect_snapshot(testout_RR$inferential$summary)
  expect_snapshot(testout_RR$inferential$fit)

  expect_snapshot(testout_RD$descriptive$summary)
  expect_snapshot(testout_RD$inferential$summary)
  expect_snapshot(testout_RD$inferential$fit)

  # Compare bootstrap outputs
  expect_snapshot(print(testout_boot_OR$descriptive$summary, digits = 5))
  expect_snapshot(print(testout_boot_OR$inferential$summary, digits = 5))
  expect_snapshot(print(testout_boot_OR$inferential$fit, digits = 5))
})
