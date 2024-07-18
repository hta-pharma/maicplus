test_that("maic_anchored works for TTE using robust SE", {
  library(flexsurv)
  ### IPD
  set.seed(1234)
  # Read in relevant ADaM data and rename variables of interest
  adsl <- read.csv(system.file("extdata", "adsl.csv",
    package = "maicplus",
    mustWork = TRUE
  ))
  adsl$USUBJID <- paste0("xx", adsl$USUBJID)
  adsl2 <- adsl
  adsl2$USUBJID <- sample(size = nrow(adsl2), paste0("yy", adsl2$USUBJID), replace = FALSE)
  adsl2 <- adsl2[order(adsl2$USUBJID), ]
  adsl <- rbind(adsl, adsl2)

  adtte <- read.csv(system.file("extdata", "adtte.csv",
    package = "maicplus",
    mustWork = TRUE
  ))
  adtte$TIME <- adtte$AVAL
  adtte$EVENT <- adtte$EVNT
  adtte$USUBJID <- paste0("xx", adtte$USUBJID)

  adtte2 <- adtte
  adtte2$ARM <- "C"
  adtte2$TIME <- adtte2$TIME * runif(nrow(adtte2), 0.15, 0.3)
  fit_C <- flexsurv::flexsurvspline(formula = Surv(TIME, EVENT) ~ 1, data = adtte2, k = 3)
  tmp <- simulate(fit_C, nsim = 1, seed = 1234, newdata = adtte2, censtime = max(adtte$TIME))
  adtte2$TIME <- tmp$time_1
  adtte2$EVENT <- tmp$event_1
  adtte2$USUBJID <- paste0("yy", adtte2$USUBJID)
  adtte <- rbind(adtte, adtte2)

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
  pseudo_ipd2 <- adtte2[, c("TIME", "EVENT", "ARM")]
  names(pseudo_ipd2) <- c("Time", "Event", "ARM")
  tmp <- simulate(fit_C, nsim = 1, seed = 4321, newdata = adtte2, censtime = max(pseudo_ipd$Time))
  pseudo_ipd2$Time <- tmp$time_1
  pseudo_ipd2$Event <- tmp$event_1
  pseudo_ipd <- rbind(pseudo_ipd, pseudo_ipd2)

  #### prepare data
  target_pop <- process_agd(target_pop)
  adsl <- dummize_ipd(adsl, dummize_cols = c("SEX"), dummize_ref_level = c("Female"))
  use_adsl <- center_ipd(ipd = adsl, agd = target_pop)

  #### derive weights
  cols <- c(
    "AGE_CENTERED", "AGE_MEDIAN_CENTERED", "AGE_SQUARED_CENTERED",
    "SEX_MALE_CENTERED", "ECOG0_CENTERED", "SMOKE_CENTERED"
  )
  # cols <- grep("_CENTERED$", names(use_adsl))
  match_res <- estimate_weights(
    data = use_adsl,
    centered_colnames = cols,
    start_val = 0,
    method = "BFGS"
  )

  # inferential result
  result <- maic_anchored(
    weights_object = match_res,
    ipd = adtte,
    trt_var_ipd = "ARM",
    pseudo_ipd = pseudo_ipd,
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
      2.56479654878613, 8.70968971110584, 2.6906650526407, 10.5753013034989,
      2.45527171390661, 4.30355056953339
    )
  )
  expect_equal(
    result$inferential$report_median_surv$`se(rmean)`,
    c(
      0.113669935856185, 0.35514766015862, 0.207503727490195, 0.573259024729393,
      0.0984888793057228, 0.336726020204787
    )
  )
  expect_equal(
    result$inferential$report_overall_robustCI$`median[95% CI]`,
    c("7.6[6.3;10.3]", "1.8[1.6; 2.0]", "12.2[10.2; NA]", " 1.8[ 1.5;2.4]", "2.7[2.3;3.3]", "1.9[1.7;2.0]", "--")
  )
  expect_equal(
    result$inferential$report_overall_robustCI$`HR[95% CI]`,
    c("0.22[0.19;0.26]", "", "0.16[0.11;0.24]", "", "0.57[0.48;0.68]", "", "0.29 [0.19; 0.44]")
  )
})


test_that("maic_anchored works for TTE using bootstrap SE", {
  # anchored example using maic_anchored for tte
  library(flexsurv)
  ### IPD
  set.seed(1234)
  # Read in relevant ADaM data and rename variables of interest
  adsl <- read.csv(system.file("extdata", "adsl.csv",
    package = "maicplus",
    mustWork = TRUE
  ))
  adsl$USUBJID <- paste0("xx", adsl$USUBJID)
  adsl2 <- adsl
  adsl2$USUBJID <- sample(size = nrow(adsl2), paste0("yy", adsl2$USUBJID), replace = FALSE)
  adsl2 <- adsl2[order(adsl2$USUBJID), ]
  adsl <- rbind(adsl, adsl2)

  adtte <- read.csv(system.file("extdata", "adtte.csv",
    package = "maicplus",
    mustWork = TRUE
  ))
  adtte$TIME <- adtte$AVAL
  adtte$EVENT <- adtte$EVNT
  adtte$USUBJID <- paste0("xx", adtte$USUBJID)

  adtte2 <- adtte
  adtte2$ARM <- "C"
  adtte2$TIME <- adtte2$TIME * runif(nrow(adtte2), 0.15, 0.3)
  fit_C <- flexsurv::flexsurvspline(formula = Surv(TIME, EVENT) ~ 1, data = adtte2, k = 3)
  tmp <- simulate(fit_C, nsim = 1, seed = 1234, newdata = adtte2, censtime = max(adtte$TIME))
  adtte2$TIME <- tmp$time_1
  adtte2$EVENT <- tmp$event_1
  adtte2$USUBJID <- paste0("yy", adtte2$USUBJID)
  adtte <- rbind(adtte, adtte2)

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
  pseudo_ipd2 <- adtte2[, c("TIME", "EVENT", "ARM")]
  names(pseudo_ipd2) <- c("Time", "Event", "ARM")
  tmp <- simulate(fit_C, nsim = 1, seed = 4321, newdata = adtte2, censtime = max(pseudo_ipd$Time))
  pseudo_ipd2$Time <- tmp$time_1
  pseudo_ipd2$Event <- tmp$event_1
  pseudo_ipd <- rbind(pseudo_ipd, pseudo_ipd2)

  #### prepare data
  target_pop <- process_agd(target_pop)
  adsl <- dummize_ipd(adsl, dummize_cols = c("SEX"), dummize_ref_level = c("Female"))
  use_adsl <- center_ipd(ipd = adsl, agd = target_pop)

  #### derive weights
  cols <- c(
    "AGE_CENTERED", "AGE_MEDIAN_CENTERED", "AGE_SQUARED_CENTERED",
    "SEX_MALE_CENTERED", "ECOG0_CENTERED", "SMOKE_CENTERED"
  )
  # cols <- grep("_CENTERED$", names(use_adsl))
  match_res_boot <- estimate_weights(
    data = use_adsl,
    centered_colnames = cols,
    start_val = 0,
    method = "BFGS",
    n_boot_iteration = 5,
    set_seed_boot = 1234
  )

  result <- maic_anchored(
    weights_object = match_res_boot,
    ipd = adtte,
    trt_var_ipd = "ARM",
    pseudo_ipd = pseudo_ipd,
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
    c("7.6[6.3;10.3]", "1.8[1.6; 2.0]", "12.2[10.2; NA]", " 1.8[ 1.5;2.4]", "2.7[2.3;3.3]", "1.9[1.7;2.0]", "--")
  )
  expect_equal(
    result$inferential$report_overall_bootCI$`HR[95% CI]`,
    c("0.22[0.19;0.26]", "", "0.16[0.11;0.24]", "", "0.57[0.48;0.68]", "", "0.29 [0.26; 0.44]")
  )

  t_matrix_expected <- matrix(
    c(
      -1.24294307251997, -1.33792890159168, -1.57070441400606, -1.52456404275395, -1.40942511608976,
      0.0466174371467349, 0.0408304991386674, 0.0405336234710769, 0.038553494493777, 0.0447691379665097,
      0.215910715682976, 0.202065581281591, 0.201329638829152, 0.196350437977044, 0.211587187623707,
      -1.80190829720242, -1.89689412627413, -2.12966963868851, -2.0835292674364, -1.96839034077221,
      0.197156067436005, 0.181888913677451, 0.181070984012272, 0.175518011252044, 0.192411579034645,
      0.0388705149268306, 0.0330835769187631, 0.0327867012511726, 0.0308065722738727, 0.0370222157466054
    ),
    byrow = FALSE,
    ncol = 6
  )
  expect_equal(result$inferential$boot_est$t, t_matrix_expected)
})


test_that("maic_anchored for binary case gives the expected result", {
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
  result <- maic_anchored(
    weights_object = weighted_twt,
    ipd = adrs_twt,
    trt_var_ipd = "ARM",
    pseudo_ipd = pseudo_adrs,
    trt_var_agd = "ARM",
    trt_ipd = "A",
    trt_agd = "B",
    trt_common = "C",
    endpoint_name = "Binary Event",
    endpoint_type = "binary",
    eff_measure = "OR"
  )

  expect_equal(
    result$inferential$report_overall_robustCI$`OR[95% CI]`,
    c("1.70[1.28;2.26]", "", "0.63[0.34;1.16]", "", "2.33[1.75;3.12]", "", "0.63 [0.34; 1.16]")
  )
  expect_equal(
    result$inferential$report_overall_robustCI$`n.events(%)`,
    c("390(78.0)", "338(67.6)", "128.8(25.8)", "115.1(23.0)", "280(58.3)", "120(37.5)", "--")
  )
  expect_equal(
    result$inferential$report_overall_robustCI$N,
    c("500", "500", "500", "500", "480", "320", "--")
  )
  expect_equal(
    result$inferential$report_overall_bootCI$`OR[95% CI]`,
    c("1.70[1.28;2.26]", "", "0.63[0.27;0.68]", "", "2.33[1.75;3.12]", "", "0.63 [0.27; 0.68]")
  )
})
