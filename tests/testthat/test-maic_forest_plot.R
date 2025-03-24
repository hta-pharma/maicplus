test_that("maic_forest_plot works for TTE", {
  data(centered_ipd_twt)
  data(adtte_twt)
  data(pseudo_ipd_twt)

  centered_colnames <- c("AGE", "AGE_SQUARED", "SEX_MALE", "ECOG0", "SMOKE", "N_PR_THER_MEDIAN")
  centered_colnames <- paste0(centered_colnames, "_CENTERED")

  #### derive weights
  weighted_data <- estimate_weights(
    data = centered_ipd_twt,
    centered_colnames = centered_colnames
  )

  # inferential result
  result <- maic_anchored(
    weights_object = weighted_data,
    ipd = adtte_twt,
    pseudo_ipd = pseudo_ipd_twt,
    trt_ipd = "A",
    trt_agd = "B",
    trt_common = "C",
    normalize_weight = FALSE,
    endpoint_name = "Overall Survival",
    endpoint_type = "tte",
    eff_measure = "HR",
    time_scale = "month",
    km_conf_type = "log-log"
  )

  make_plot <- function() {maic_forest_plot(
    result,
    xlim = c(0, 3.5),
    reference_line = 1)
  }
  expect_error(make_plot(), NA)

  skip_if_not_installed("vdiffr")
  vdiffr::expect_doppelganger(
    title = "maic_forest_plot_tte",
    fig = make_plot()
  )
  })
