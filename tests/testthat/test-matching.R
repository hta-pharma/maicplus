test_that("ess_footnote_text works", {
  expect_equal(
    ess_footnote_text(width = 100),
    paste0(
      "An ESS reduction up to ~60% is not unexpected based on the 2021 survey of NICE's technology\n",
      "appraisals (https://onlinelibrary.wiley.com/doi/full/10.1002/jrsm.1511), whereas a reduction of\n",
      ">75% is less common and it may be considered suboptimal."
    )
  )
})

test_that("calculate_weights_legend works", {
  load(system.file("extdata", "weighted_data.rda", package = "maicplus", mustWork = TRUE))
  result <- calculate_weights_legend(weighted_data)
  expect_equal(
    result,
    expected = list(
      ess = 166.37,
      ess_reduction = 66.73,
      wt_median = 0.0594,
      wt_scaled_median = 0.1486,
      nr_na = 0L
    )
  )
})

test_that("optimise_weights works as expected", {
  set.seed(123)
  object <- matrix(c(age_centered = rnorm(20, sd = 2), biomarker_centered = rnorm(20, sd = 1)), nrow = 20)
  expect_output(
    result <- optimise_weights(object, par = c(0, 0), trace = 2),
    "converged"
  )
  expect_equal(result$alpha, c(-0.0770222925868337, 0.0682633331364419))
  expect_equal(
    result$wt[, 1],
    c(
      1.01353582399708, 1.02079108647525, 0.733337637418428, 0.941182721058724,
      0.9393347760529, 0.684315793280964, 0.986283431666655, 1.22794936276493,
      1.0285118397984, 1.16677251045601, 0.852611485594579, 0.927214883798323,
      0.999367789272156, 1.04382749980235, 1.15224035280605, 0.795920944660199,
      0.961867179319782, 1.34813266573421, 0.879038915192837, 1.04797395929194
    )
  )
})

test_that("estimate_weights fails as expected with invalid input", {
  load(system.file("extdata", "ipd.rda", package = "maicplus", mustWork = TRUE))
  load(system.file("extdata", "agd.rda", package = "maicplus", mustWork = TRUE))
  ipd_centered <- center_ipd(ipd = ipd, agd = agd)
  centered_colnames <- paste0(c("AGE", "AGE_SQUARED", "SEX_MALE", "ECOG0", "SMOKE", "N_PR_THER_MEDIAN"), "_CENTERED")
  expect_error(
    estimate_weights(data = as.matrix(ipd_centered), centered_colnames = centered_colnames),
    "not a data.frame"
  )

  expect_error(
    estimate_weights(data = ipd_centered, centered_colnames = "hba1c_centered", trace = 2),
    "specified centered_colnames"
  )

  ipd_centered$AGE_MEDIAN_CENTERED_FCT <- as.factor(ipd_centered$AGE_MEDIAN_CENTERED)
  expect_error(
    estimate_weights(data = ipd_centered, centered_colnames = c(centered_colnames, "AGE_MEDIAN_CENTERED_FCT")),
    "not numeric"
  )
})



test_that("estimate_weights works as expected", {
  load(system.file("extdata", "ipd.rda", package = "maicplus", mustWork = TRUE))
  load(system.file("extdata", "agd.rda", package = "maicplus", mustWork = TRUE))
  ipd_centered <- center_ipd(ipd = ipd, agd = agd)
  centered_colnames <- paste0(c("AGE", "AGE_SQUARED", "SEX_MALE", "ECOG0", "SMOKE", "N_PR_THER_MEDIAN"), "_CENTERED")
  expect_output(
    result <- estimate_weights(data = ipd_centered, centered_colnames = centered_colnames, trace = 2),
    "converged"
  )

  expect_s3_class(result, "maicplus_estimate_weights")
  expect_equal(sum(result$data$weights), 199.8422368)
  expect_equal(sum(result$data$scaled_weights), 500)
  expect_equal(result$ess, 166.3675302)
  expect_null(result$boot)
})

test_that("estimate_weights works as expected with bootstrapping", {
  load(system.file("extdata", "ipd.rda", package = "maicplus", mustWork = TRUE))
  load(system.file("extdata", "agd.rda", package = "maicplus", mustWork = TRUE))
  ipd_centered <- center_ipd(ipd = ipd, agd = agd)
  centered_colnames <- paste0(c("AGE", "AGE_SQUARED", "SEX_MALE", "ECOG0", "SMOKE", "N_PR_THER_MEDIAN"), "_CENTERED")
  expect_output(
    result <- estimate_weights(
      data = ipd_centered,
      centered_colnames = centered_colnames,
      n_boot_iteration = 3,
      set_seed = 999,
      trace = 2
    ),
    "converged"
  )

  expect_s3_class(result, "maicplus_estimate_weights")
  expect_equal(sum(result$data$weights), 199.8422368)
  expect_equal(sum(result$data$scaled_weights), 500)
  expect_equal(result$ess, 166.3675302)
  expect_equal(dim(result$boot), c(500, 2, 3))

  expected_matrix <- structure(
    c(
      c(411, 71, 321, 150, 440),
      c(8.71745219703408e-05, 1.93529165335151e-08, 0.886392116942025, 1.42604077815808e-13, 0.802986096194834)
    ),
    dim = c(5L, 2L),
    dimnames = list("sampled_patient" = NULL, c("rowid", "weight"))
  )
  expect_equal(result$boot[1:5, , 1], expected_matrix)
})


test_that("estimate_weights works as expected with alternative bootstrap strata", {
  load(system.file("extdata", "ipd.rda", package = "maicplus", mustWork = TRUE))
  load(system.file("extdata", "agd.rda", package = "maicplus", mustWork = TRUE))
  ipd_centered <- center_ipd(ipd = ipd, agd = agd)
  centered_colnames <- paste0(c("AGE", "AGE_SQUARED", "ECOG0", "SMOKE", "N_PR_THER_MEDIAN"), "_CENTERED")
  expect_output(
    result <- estimate_weights(
      data = ipd_centered,
      centered_colnames = centered_colnames,
      n_boot_iteration = 3,
      set_seed = 999,
      trace = 2,
      boot_strata = c("ARM", "SEX")
    ),
    "converged"
  )

  expect_s3_class(result, "maicplus_estimate_weights")
  expect_equal(sum(result$data$weights), 206.83843133)

  expect_error(
    result <- estimate_weights(
      data = ipd_centered,
      centered_colnames = centered_colnames,
      n_boot_iteration = 3,
      set_seed = 999,
      trace = 2,
      boot_strata = "FISH"
    ),
    "boot_strata are not in data"
  )
})


test_that("estimate_weights works when the input is a tibble", {
  skip_if_not_installed("tibble")
  load(system.file("extdata", "ipd.rda", package = "maicplus", mustWork = TRUE))
  load(system.file("extdata", "agd.rda", package = "maicplus", mustWork = TRUE))
  ipd_centered <- center_ipd(ipd = ipd, agd = agd)
  centered_colnames <- paste0(c("AGE", "AGE_SQUARED", "SEX_MALE", "ECOG0", "SMOKE", "N_PR_THER_MEDIAN"), "_CENTERED")
  ipd_centered <- tibble::as_tibble(ipd_centered)
  expect_output(
    result <- estimate_weights(data = ipd_centered, centered_colnames = centered_colnames, trace = 2),
    "converged"
  )

  expect_s3_class(result, "maicplus_estimate_weights")
  expect_equal(sum(result$data$weights), 199.8422368)
  expect_equal(sum(result$data$scaled_weights), 500)
  expect_equal(result$ess, 166.3675302)
  expect_null(result$boot)
})


test_that("estimate_weights prints errors about convergence", {
  load(system.file("extdata", "ipd.rda", package = "maicplus", mustWork = TRUE))
  load(system.file("extdata", "agd.rda", package = "maicplus", mustWork = TRUE))
  ipd_centered <- center_ipd(ipd = ipd, agd = agd)
  centered_colnames <- paste0(c("AGE", "AGE_SQUARED", "SEX_MALE", "ECOG0", "SMOKE", "N_PR_THER_MEDIAN"), "_CENTERED")
  expect_warning(
    estimate_weights(data = ipd_centered, centered_colnames = centered_colnames, trace = 0, maxit = 3),
    "did not converge"
  )

  expect_output(
    estimate_weights(data = ipd_centered, centered_colnames = centered_colnames, trace = 0, maxit = 300),
    NA
  )
})

test_that("plot_weights_base works as expected", {
  vdiffr::expect_doppelganger(
    title = "plot_weights_base scaled_TRUE",
    plot_weights_base(
      weighted_sat,
      bin_col = "#6ECFFF",
      vline_col = "#0000E8",
      main_title = c("Scaled Individual Weights"),
      scaled_weights = TRUE
    )
  )

  vdiffr::expect_doppelganger(
    title = "plot_weights_base scaled_FALSE",
    plot_weights_base(
      weighted_sat,
      bin_col = "#6ECFFF",
      vline_col = "#0000E8",
      main_title = c("Unscaled Individual Weights"),
      scaled_weights = FALSE
    )
  )
})

test_that("plot_weights_ggplot works as expected", {
  vdiffr::expect_doppelganger(
    title = "plot_weights_ggplot",
    plot_weights_ggplot(
      weighted_sat,
      bin_col = "#6ECFFF",
      vline_col = "#0000E8",
      main_title = c("Scaled Individual Weights", "Unscaled Individual Weights"),
      bins = 10
    )
  )
})

test_that("default plot works as expected", {
  vdiffr::expect_doppelganger(
    title = "default_weights_plot",
    plot(weighted_twt)
  )
  vdiffr::expect_doppelganger(
    title = "default_weights_ggplot",
    plot(weighted_twt, ggplot = TRUE)
  )
})

test_that("check_weights works as expected", {
  result <- check_weights(weighted_sat, process_agd(agd))
  checkmate::expect_class(result, "maicplus_check_weights")
  checkmate::expect_data_frame(result, nrows = 6, ncols = 6)
  checkmate::expect_names(
    colnames(result),
    identical.to = c(
      "covariate", "match_stat", "internal_trial", "internal_trial_after_weighted",
      "external_trial", "sum_centered_IPD_with_weights"
    ),
    what = "colnames"
  )
  expect_equal(result[, "covariate"], c("AGE", "AGE", "AGE", "SEX_MALE", "ECOG0", "SMOKE"))
  expect_equal(result[, "match_stat"], c("Mean", "Median", "SD", "Prop", "Prop", "Prop"))
  expect_equal(result[, "sum_centered_IPD_with_weights"], c(0, 0, -0.0045, 0, 0, 0))
})
