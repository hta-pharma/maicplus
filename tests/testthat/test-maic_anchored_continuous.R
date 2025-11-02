test_that("maic_anchored works with continuous MD endpoint", {
  skip_if_not_installed("maicplus")
  
  data("centered_ipd_twt", package = "maicplus")
  data("adlb_twt", package = "maicplus")
  
  centered_colnames <- paste0(c("AGE", "SEX_MALE", "ECOG0"), "_CENTERED")
  
  weighted_data <- estimate_weights(
    data = centered_ipd_twt,
    centered_colnames = centered_colnames
  )
  
  pseudo_adlb <- get_pseudo_ipd_continuous(
    data.frame(ARM = c("B", "C"), MEAN = c(-12, -9), SD = c(9, 9.5), N = c(120, 150))
  )
  
  result <- maic_anchored(
    weights_object = weighted_data,
    ipd = adlb_twt,
    pseudo_ipd = pseudo_adlb,
    trt_ipd = "A",
    trt_agd = "B",
    trt_common = "C",
    endpoint_type = "continuous",
    eff_measure = "MD"
  )
  
  expect_type(result, "list")
  expect_true("descriptive" %in% names(result))
  expect_true("inferential" %in% names(result))
  expect_s3_class(result$inferential$summary, "data.frame")
})

test_that("maic_anchored works with continuous SMD endpoint", {
  skip_if_not_installed("maicplus")
  
  data("centered_ipd_twt", package = "maicplus")
  data("adlb_twt", package = "maicplus")
  
  centered_colnames <- paste0(c("AGE", "SEX_MALE"), "_CENTERED")
  
  weighted_data <- estimate_weights(
    data = centered_ipd_twt,
    centered_colnames = centered_colnames
  )
  
  pseudo_adlb <- get_pseudo_ipd_continuous(
    data.frame(ARM = c("B", "C"), MEAN = c(-12, -9), SD = c(9, 9.5), N = c(120, 150))
  )
  
  result <- maic_anchored(
    weights_object = weighted_data,
    ipd = adlb_twt,
    pseudo_ipd = pseudo_adlb,
    trt_ipd = "A",
    trt_agd = "B",
    trt_common = "C",
    endpoint_type = "continuous",
    eff_measure = "SMD",
    smd_type = "cohen"
  )
  
  expect_type(result, "list")
  expect_true(all(c("descriptive", "inferential") %in% names(result)))
})

