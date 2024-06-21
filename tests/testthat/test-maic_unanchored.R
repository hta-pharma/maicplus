test_that("test binary case", {
  # load in prognostic IPD data and AgD
  load(system.file("extdata", "ipd.rda", package = "maicplus", mustWork = TRUE))
  load(system.file("extdata", "agd.rda", package = "maicplus", mustWork = TRUE))
  ipd_centered <- center_ipd(ipd = ipd, agd = agd)

  # estimate weights
  centered_colnames <- c("AGE", "AGE_SQUARED", "SEX_MALE", "ECOG0", "SMOKE", "N_PR_THER_MEDIAN")
  centered_colnames <- paste0(centered_colnames, "_CENTERED")

  weighted_data <- estimate_weights(data = ipd_centered, centered_colnames = centered_colnames)
  weighted_data2 <- estimate_weights(
    data = ipd_centered, centered_colnames = centered_colnames, n_boot_iteration = 20,
    set_seed_boot = 1234
  )

  # get dummy binary IPD
  adrs <- read.csv(system.file("extdata", "adrs.csv", package = "maicplus", mustWork = TRUE))
  adrs$RESPONSE <- adrs$AVAL

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
      ipd = adrs,
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
      ipd = adrs,
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
  # anchored example using maic_anchored for tte
  # library(flexsurv)

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
  cols <- c(
    "AGE_CENTERED", "AGE_MEDIAN_CENTERED", "AGE_SQUARED_CENTERED",
    "SEX_MALE_CENTERED", "ECOG0_CENTERED", "SMOKE_CENTERED"
  )
  # cols <-  grep("_CENTERED$", names(use_adsl))
  match_res <- estimate_weights(
    data = use_adsl,
    centered_colnames = cols,
    start_val = 0,
    method = "BFGS"
  )

  match_res_boot <- estimate_weights(
    data = use_adsl,
    centered_colnames = cols,
    start_val = 0,
    method = "BFGS",
    n_boot_iteration = 500,
    set_seed_boot = 1234
  )

  # inferential result
  testout <- maic_unanchored(
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

  testout2 <- maic_unanchored(
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
