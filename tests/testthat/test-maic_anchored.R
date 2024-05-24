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
  match_res <- estimate_weights(
    data = use_adsl,
    centered_colnames = grep("_CENTERED$", names(use_adsl)),
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
    c(
      "7.6[6.3;10.3]", "1.8[1.6; 2.0]", "12.2[10.2; NA]", " 1.8[ 1.5;2.4]",
      "2.7[2.3;3.3]", "1.9[1.7;2.0]", "-"
    )
  )
  expect_equal(
    result$inferential$report_overall_robustCI$`HR[95% CI]`,
    c(
      "0.22[0.19;0.26]", "", "0.16[0.11;0.24]", "", "0.57[0.48;0.68]",
      "", "0.29 [0.19; 0.44]"
    )
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
  match_res_boot <- estimate_weights(
    data = use_adsl,
    centered_colnames = grep("_CENTERED$", names(use_adsl)),
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
    c(
      "7.6[6.3;10.3]", "1.8[1.6; 2.0]", "12.2[10.2; NA]", " 1.8[ 1.5;2.4]",
      "2.7[2.3;3.3]", "1.9[1.7;2.0]", "-"
    )
  )
  expect_equal(
    result$inferential$report_overall_bootCI$`HR[95% CI]`,
    c(
      "0.22[0.19;0.26]", "", "0.16[0.11;0.24]", "", "0.57[0.48;0.68]",
      "", "0.29 [0.20; 0.38]"
    )
  )

  expect_equal(
    result$inferential$boot_est$t,
    matrix(
      c(
        -0.473643217943453, -0.0280612969274431, -0.292544997523232, -0.257912426911704, -0.332304813173113,
        0.0243022844964305, 0.0180627341532462, 0.0160751359632564, 0.0182590931935275, 0.0215193429847375
      ),
      nrow = 5
    )
  )
})
