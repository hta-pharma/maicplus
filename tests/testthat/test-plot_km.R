test_that("kmplot works by trial", {
  make_plot <- function() {
    kmplot(
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
      time_grid = seq(0, 20, by = 2),
      use_colors = NULL,
      use_line_types = NULL,
      use_pch_cex = 0.65,
      use_pch_alpha = 100
    )
  }

  expect_no_error(make_plot())
  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger(
    title = "kmplot by_trial",
    fig = make_plot()
  )
})

test_that("kmplot works all", {
  make_plot <- function() {
    kmplot(
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
      time_scale = "year",
      time_grid = seq(0, 1.5, by = 0.25),
      use_colors = NULL,
      use_line_types = c(6, 3, 4, 5),
      use_pch_cex = 0.65,
      use_pch_alpha = 100
    )
  }

  expect_no_error(make_plot())
  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger(
    title = "kmplot_all",
    fig = make_plot()
  )
})


test_that("kmplot normalize weights works", {
  make_plot <- function() {
    kmplot(
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
      km_layout = "all",
      time_scale = "year",
      time_grid = seq(0, 1.5, by = 0.25),
      use_colors = NULL,
      use_line_types = c(6, 3, 4, 5),
      use_pch_cex = 0.65,
      use_pch_alpha = 100
    )
  }

  expect_no_error(make_plot())
  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger(
    title = "kmplot_normalize_all",
    fig = make_plot()
  )
})

test_that("kmplot works by_arm", {
  make_plot <- function() {
    kmplot(
      weights_object = weighted_twt,
      tte_ipd = adtte_twt,
      tte_pseudo_ipd = pseudo_ipd_twt,
      trt_ipd = "A",
      trt_agd = "B",
      trt_common = "C",
      trt_var_ipd = "ARM",
      trt_var_agd = "ARM",
      endpoint_name = "PFS",
      km_conf_type = "log-log",
      km_layout = "by_arm",
      use_colors = c("red", "yellow", "orange"),
      use_line_types = NULL,
      use_pch_cex = 0.75,
      use_pch_alpha = 85
    )
  }

  expect_no_error(make_plot())
  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger(
    title = "kmplot by_arm",
    fig = make_plot()
  )
})



test_that("ph_diagplot works for anchored data", {
  expect_no_error(
    ph_diagplot(
      weights_object = weighted_twt,
      tte_ipd = adtte_twt,
      tte_pseudo_ipd = pseudo_ipd_twt,
      trt_ipd = "A",
      trt_agd = "B",
      trt_common = "C",
      trt_var_ipd = "ARM",
      trt_var_agd = "ARM",
      endpoint_name = "Time to Event Endpoint",
      time_scale = "months",
      zph_transform = "log",
      zph_log_hazard = TRUE
    )
  )

  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger(
    title = "ph_diagplot anchored",
    fig = function() {
      ph_diagplot(
        weights_object = weighted_twt,
        tte_ipd = adtte_twt,
        tte_pseudo_ipd = pseudo_ipd_twt,
        trt_ipd = "A",
        trt_agd = "B",
        trt_common = "C",
        trt_var_ipd = "ARM",
        trt_var_agd = "ARM",
        endpoint_name = "Time to Event Endpoint",
        time_scale = "years",
        zph_transform = "log",
        zph_log_hazard = TRUE
      )
    }
  )
})

test_that("ph_diagplot works for unanchored data", {
  expect_no_error(
    ph_diagplot(
      weights_object = weighted_sat,
      tte_ipd = adtte_sat,
      tte_pseudo_ipd = pseudo_ipd_sat,
      trt_ipd = "A",
      trt_agd = "B",
      trt_var_ipd = "ARM",
      trt_var_agd = "ARM",
      endpoint_name = "Time to Event Endpoint",
      time_scale = "days",
      zph_transform = "log",
      zph_log_hazard = FALSE
    )
  )

  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger(
    title = "ph_diagplot unanchored",
    fig = function() {
      ph_diagplot(
        weights_object = weighted_twt,
        tte_ipd = adtte_twt,
        tte_pseudo_ipd = pseudo_ipd_twt,
        trt_ipd = "A",
        trt_agd = "B",
        trt_common = "C",
        trt_var_ipd = "ARM",
        trt_var_agd = "ARM",
        endpoint_name = "Time to Event Endpoint",
        time_scale = "weeks",
        zph_transform = "log",
        zph_log_hazard = TRUE
      )
    }
  )
})

test_that("ph_diagplot_lch works without error", {
  original_par <- par("mar", "bty", "tcl", "mgp")
  data(adtte_sat)
  data(pseudo_ipd_sat)
  combined_data <- rbind(adtte_sat[, c("TIME", "EVENT", "ARM")], pseudo_ipd_sat)
  kmobj <- survfit(Surv(TIME, EVENT) ~ ARM, combined_data, conf.type = "log-log")
  expect_no_error(
    ph_diagplot_lch(
      kmobj,
      time_scale = "month", log_time = TRUE,
      endpoint_name = "OS", subtitle = "(Before Matching)"
    )
  )

  expect_equal(original_par, par(names(original_par)))

  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger(
    title = "ph_diagplot_lch",
    fig = function() {
      ph_diagplot_lch(
        kmobj,
        time_scale = "month", log_time = TRUE,
        endpoint_name = "OS", subtitle = "(Before Matching)"
      )
    }
  )
})

test_that("ph_diagplot_schoenfeld works without error", {
  original_par <- par("bty", "mar", "tcl", "mgp")
  data(adtte_sat)
  data(pseudo_ipd_sat)
  combined_data <- rbind(adtte_sat[, c("TIME", "EVENT", "ARM")], pseudo_ipd_sat)
  unweighted_cox <- coxph(Surv(TIME, EVENT == 1) ~ ARM, data = combined_data)
  expect_no_error(
    ph_diagplot_schoenfeld(
      unweighted_cox,
      time_scale = "month", log_time = TRUE,
      endpoint_name = "OS", subtitle = "(Before Matching)"
    )
  )
  expect_equal(original_par, par(names(original_par)))

  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger(
    title = "ph_diagplot_schoenfeld",
    fig = function() {
      ph_diagplot_schoenfeld(
        unweighted_cox,
        time_scale = "month", log_time = TRUE,
        endpoint_name = "OS", subtitle = "(Before Matching)"
      )
    }
  )
})
