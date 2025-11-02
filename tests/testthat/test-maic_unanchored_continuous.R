test_that("maic_unanchored works with continuous MD endpoint", {
  skip_if_not_installed("maicplus")
  
  data("centered_ipd_sat", package = "maicplus")
  data("adlb_sat", package = "maicplus")
  
  centered_colnames <- paste0(c("AGE", "SEX_MALE", "ECOG0"), "_CENTERED")
  
  weighted_data <- estimate_weights(
    data = centered_ipd_sat,
    centered_colnames = centered_colnames
  )
  
  pseudo_adlb <- get_pseudo_ipd_continuous(
    data.frame(ARM = "B", MEAN = -12, SD = 9, N = 120)
  )
  
  result <- maic_unanchored(
    weights_object = weighted_data,
    ipd = adlb_sat,
    pseudo_ipd = pseudo_adlb,
    trt_ipd = "A",
    trt_agd = "B",
    endpoint_type = "continuous",
    eff_measure = "MD"
  )
  
  expect_type(result, "list")
  expect_true("descriptive" %in% names(result))
  expect_true("inferential" %in% names(result))
  expect_s3_class(result$inferential$summary, "data.frame")
})

test_that("maic_unanchored works with continuous SMD endpoint", {
  skip_if_not_installed("maicplus")
  
  data("centered_ipd_sat", package = "maicplus")
  data("adlb_sat", package = "maicplus")
  
  centered_colnames <- paste0(c("AGE", "SEX_MALE"), "_CENTERED")
  
  weighted_data <- estimate_weights(
    data = centered_ipd_sat,
    centered_colnames = centered_colnames
  )
  
  pseudo_adlb <- get_pseudo_ipd_continuous(
    data.frame(ARM = "B", MEAN = -12, SD = 9, N = 120)
  )
  
  result <- maic_unanchored(
    weights_object = weighted_data,
    ipd = adlb_sat,
    pseudo_ipd = pseudo_adlb,
    trt_ipd = "A",
    trt_agd = "B",
    endpoint_type = "continuous",
    eff_measure = "SMD",
    smd_type = "hedges"
  )
  
  expect_type(result, "list")
  expect_true(all(c("descriptive", "inferential") %in% names(result)))
})

test_that("maic_unanchored validates RoM positive values", {
  skip_if_not_installed("maicplus")
  
  data("centered_ipd_sat", package = "maicplus")
  
  centered_colnames <- paste0(c("AGE", "SEX_MALE"), "_CENTERED")
  
  weighted_data <- estimate_weights(
    data = centered_ipd_sat,
    centered_colnames = centered_colnames
  )
  
  # Create negative values using SAME USUBJIDs as weighted_data
  adlb_negative <- data.frame(
    USUBJID = weighted_data$data$USUBJID,
    ARM = "A",
    AVAL = rnorm(nrow(weighted_data$data), -10, 2)
  )
  
  pseudo_adlb <- get_pseudo_ipd_continuous(
    data.frame(ARM = "B", MEAN = 10, SD = 2, N = 100)
  )
  
  expect_error(
    maic_unanchored(
      weights_object = weighted_data,
      ipd = adlb_negative,
      pseudo_ipd = pseudo_adlb,
      trt_ipd = "A",
      trt_agd = "B",
      endpoint_type = "continuous",
      eff_measure = "RoM"
    ),
    "RoM.*positive"
  )
})

