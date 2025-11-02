test_that("get_pseudo_ipd_continuous works correctly", {
  # Test basic functionality
  agd <- data.frame(
    ARM = "B",
    MEAN = 50,
    SD = 10,
    N = 100
  )
  
  result <- get_pseudo_ipd_continuous(agd)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 100)
  expect_true(all(c("USUBJID", "ARM", "AVAL") %in% names(result)))
  expect_equal(unique(result$ARM), "B")
})

test_that("get_pseudo_ipd_continuous handles multiple arms", {
  agd <- data.frame(
    ARM = c("B", "C"),
    MEAN = c(50, 55),
    SD = c(10, 12),
    N = c(100, 120)
  )
  
  result <- get_pseudo_ipd_continuous(agd)
  
  expect_equal(nrow(result), 220)
  expect_setequal(unique(result$ARM), c("B", "C"))
})

test_that("get_pseudo_ipd_continuous validates inputs", {
  expect_error(get_pseudo_ipd_continuous("not a dataframe"))
  expect_error(get_pseudo_ipd_continuous(data.frame(A = 1)))
  expect_error(get_pseudo_ipd_continuous(data.frame(ARM = "A", MEAN = "x", SD = 1, N = 10)))
  expect_error(get_pseudo_ipd_continuous(data.frame(ARM = "A", MEAN = 1, SD = -1, N = 10)))
  expect_error(get_pseudo_ipd_continuous(data.frame(ARM = "A", MEAN = 1, SD = 1, N = -10)))
})

test_that("lm_makeup works for unweighted models", {
  data("adlb_twt", package = "maicplus")
  # adlb_twt has arms A and C
  test_data <- adlb_twt[adlb_twt$ARM %in% c("A", "C"), ]
  test_data$ARM <- factor(test_data$ARM, levels = c("A", "C"))
  lmobj <- lm(AVAL ~ ARM, data = test_data)
  
  result <- lm_makeup(lmobj, legend = "test")
  
  expect_s3_class(result, "data.frame")
  expect_true(all(c("treatment", "type", "n", "mean", "sd") %in% names(result)))
  expect_equal(result$type[1], "test")
})

test_that("calculate_smd computes correctly for Cohen's d", {
  # Access internal function
  smd <- maicplus:::calculate_smd(mean1 = 10, mean2 = 8, sd1 = 2, sd2 = 2, n1 = 100, n2 = 100, method = "cohen")
  expect_type(smd, "double")
  expect_equal(smd, 1, tolerance = 0.01)  # (10-8)/2 = 1
})

test_that("calculate_smd methods differ appropriately", {
  smd_cohen <- maicplus:::calculate_smd(10, 8, 2, 3, 100, 100, "cohen")
  smd_hedges <- maicplus:::calculate_smd(10, 8, 2, 3, 100, 100, "hedges")
  smd_glass <- maicplus:::calculate_smd(10, 8, 2, 3, 100, 100, "glass")
  
  expect_true(smd_cohen != smd_glass)  # Glass uses only sd2
  expect_true(abs(smd_cohen - smd_hedges) < 0.1)  # Hedges is corrected Cohen
})

test_that("calculate_smd_se computes standard errors", {
  se <- maicplus:::calculate_smd_se(smd = 0.5, n1 = 100, n2 = 100, method = "cohen")
  expect_type(se, "double")
  expect_true(se > 0)
})

