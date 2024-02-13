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
    result <- optimise_weights(object, par = c(0, 0)),
    "converged"
  )
  expect_equal(result$alpha, c(-0.07702229, 0.06826333))
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


test_that("estimate_weights works as expected", {
  load(system.file("extdata", "ipd.rda", package = "maicplus", mustWork = TRUE))
  load(system.file("extdata", "agd.rda", package = "maicplus", mustWork = TRUE))
  ipd_centered <- center_ipd(ipd = ipd, agd = agd)
  centered_colnames <- paste0(c("AGE", "AGE_SQUARED", "SEX_MALE", "ECOG0", "SMOKE", "N_PR_THER_MEDIAN"), "_CENTERED")
  expect_output(
    result <- estimate_weights(data = ipd_centered, centered_colnames = centered_colnames),
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
      set_seed = 999
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
      c(411, 324, 61, 71, 361),
      c(0.0001657791105167, 0.382238891368213, 1.41902681600543e-10, 5.04595364860806e-08, 0.434732317134287)
    ),
    dim = c(5L, 2L),
    dimnames = list(c("411", "324", "61", "71", "361"), c("rowid", "weight"))
  )
  expect_equal(result$boot[1:5, , 1], expected_matrix)
})
