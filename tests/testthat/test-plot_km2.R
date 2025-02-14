test_that("kmplot2 works", {
  skip_if_not_installed("vdiffr")

  data(weighted_twt)
  data(adtte_twt)
  data(pseudo_ipd_twt)

  vdiffr::expect_doppelganger(
    title = "kmplot2_by_trial",
    fig = function() {
      suppressWarnings(
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
      )
    }
  )

  vdiffr::expect_doppelganger(
    title = "kmplot2_normalize_by_trial",
    fig = function() {
      kmplot2(
        weights_object = weighted_twt,
        tte_ipd = adtte_twt,
        tte_pseudo_ipd = pseudo_ipd_twt,
        trt_ipd = "A",
        trt_agd = "B",
        trt_common = "C",
        normalize_weights = TRUE,
        trt_var_ipd = "ARM",
        trt_var_agd = "ARM",
        endpoint_name = "Overall Survival",
        km_conf_type = "log-log",
        km_layout = "by_trial",
        time_scale = "month",
        break_x_by = 2
      )
    }
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
