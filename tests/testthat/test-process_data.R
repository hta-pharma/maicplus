test_that("process_agd works", {
  raw_agd <- maicplus::agd
  result <- process_agd(raw_agd)
  expected <- data.frame(
    STUDY = "Study_XXXX",
    ARM = "Total",
    N = 300,
    AGE_MEAN = 51,
    AGE_MEDIAN = 49,
    AGE_SD = 3.25,
    SEX_MALE_PROP = 0.49,
    ECOG0_PROP = 0.35,
    SMOKE_PROP = 0.19333333333333,
    N_PR_THER_MEDIAN = 2
  )
  expect_equal(result, expected)

  colnames(raw_agd)[2] <- "TREATMENT"
  expect_warning(
    result2 <- process_agd(raw_agd),
    "renamed"
  )
  expect_equal(colnames(result2), colnames(result))
})

test_that("dummize_ipd works as expected", {
  dat <- maicplus::adsl_twt
  result <- dummize_ipd(dat, dummize_cols = c("SEX"), dummize_ref_level = c("Male"))
  checkmate::expect_data_frame(result, nrow = 1000, ncol = 9)
  expect_equal(
    colnames(result),
    c("USUBJID", "ARM", "AGE", "SEX", "SMOKE", "ECOG0", "N_PR_THER", "SEX_MALE", "SEX_FEMALE")
  )

  result2 <- dummize_ipd(dat, dummize_cols = c("SEX", "SMOKE"), dummize_ref_level = c("Female", "1"))
  checkmate::expect_data_frame(result2, nrow = 1000, ncol = 10)

  # Note: this data already has a "SEX_MALE" column, so that's why we have 2 in this result
  expect_equal(
    colnames(result2),
    c("USUBJID", "ARM", "AGE", "SEX", "SMOKE", "ECOG0", "N_PR_THER", "SEX_MALE", "SEX_MALE", "SMOKE_0")
  )
})

test_that("complete_agd works as expected", {
  dat <- maicplus::agd[, -c(5, 7, 8)]
  dat[, "N"] <- c(NA, 200, 100)
  dat[, "AGE_MEAN"] <- c(NA, 50, 60)
  dat[, "AGE_SD"] <- c(NA, 6, 7)
  dat[, "SMOKE_COUNT"] <- c(NA, 10, 20)
  dat[, "N_PR_THER_MEDIAN"] <- c(NA, 3, 4)

  result <- complete_agd(dat)

  expected <- data.frame(
    STUDY = c("Study_XXXX", "Study_XXXX", "Study_XXXX"),
    ARM = c("Intervention", "Comparator", "total"),
    N = c(200, 100, 300),
    AGE_MEAN = c(50, 60, 53.3333333333333),
    AGE_SD = c(6, 7, 6.33908088671333),
    SMOKE_COUNT = c(10, 20, 30),
    N_PR_THER_MEDIAN = c(3, 4, 3.5)
  )

  expect_equal(result, expected)
  expect_equal(result[3, "AGE_MEAN"], (50 * 200 + 100 * 60) / 300)
  expect_equal(result[3, "AGE_SD"], sqrt(((200 - 1) * 6^2 + (100 - 1) * 7^2) / (300 - 1)))
})
