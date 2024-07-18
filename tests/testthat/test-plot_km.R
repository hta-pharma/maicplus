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
