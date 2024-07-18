test_that("ph_diagplot_lch works without error", {
  original_par <- par()
  load(system.file("extdata", "combined_data_tte.rda", package = "maicplus", mustWork = TRUE))
  kmobj <- survfit(Surv(TIME, EVENT) ~ ARM, combined_data_tte, conf.type = "log-log")
  ph_diagplot_lch(kmobj,
    time_scale = "month", log_time = TRUE,
    endpoint_name = "OS", subtitle = "(Before Matching)"
  )

  expect_equal(original_par, par())
})

test_that("ph_diagplot_schoenfeld works without error", {
  original_par <- par()
  load(system.file("extdata", "combined_data_tte.rda", package = "maicplus", mustWork = TRUE))
  unweighted_cox <- coxph(Surv(TIME, EVENT == 1) ~ ARM, data = combined_data_tte)
  expect_no_error(
    ph_diagplot_schoenfeld(unweighted_cox,
      time_scale = "month", log_time = TRUE,
      endpoint_name = "OS", subtitle = "(Before Matching)"
    )
  )
  expect_equal(original_par, par())
})

test_that("kmplot2 works", {
  skip_if_not_installed("vdiffr")

  data(weighted_twt)
  data(adtte_twt)
  data(pseudo_ipd_twt)

  # plot by trial
  expect_warning(
    vdiffr::expect_doppelganger(
      title = "kmplot2_by_trial",
      fig = function() {
        kmplot2(
          weights_object = weighted_twt,
          tte_ipd = adtte_twt,
          tte_pseudo_ipd = pseudo_ipd_twt,
          trt_ipd = "A",
          trt_agd = "B",
          trt_common = "C",
          trt_var_ipd = "ARM",
          trt_var_agd = "ARM",
          endpoint_name = "Overall Survival",
          km_conf_type = "log-log",
          km_layout = "by_trial",
          time_scale = "month",
          break_x_by = 2
        )
      }
    ),
    "select_"
  )



  # plot by arm
  vdiffr::expect_doppelganger(
    title = "kmplot2_by_am",
    fig = function() {
      kmplot2(
        weights_object = weighted_twt,
        tte_ipd = adtte_twt,
        tte_pseudo_ipd = pseudo_ipd_twt,
        trt_ipd = "A",
        trt_agd = "B",
        trt_common = "C",
        trt_var_ipd = "ARM",
        trt_var_agd = "ARM",
        endpoint_name = "Overall Survival",
        km_conf_type = "log-log",
        km_layout = "by_arm",
        time_scale = "month",
        break_x_by = 2
      )
    }
  )

  # plot all
  vdiffr::expect_doppelganger(
    title = "kmplot2_all",
    fig = function() {
      kmplot2(
        weights_object = weighted_twt,
        tte_ipd = adtte_twt,
        tte_pseudo_ipd = pseudo_ipd_twt,
        trt_ipd = "A",
        trt_agd = "B",
        trt_common = "C",
        trt_var_ipd = "ARM",
        trt_var_agd = "ARM",
        endpoint_name = "Overall Survival",
        km_conf_type = "log-log",
        km_layout = "all",
        time_scale = "month",
        break_x_by = 2,
        xlim = c(0, 20),
        show_risk_set = FALSE
      )
    }
  )
})
