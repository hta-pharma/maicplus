test_that("generate pseudo binary IPD given unstacked table as input", {
  # example of unstacked
  testdat <- data.frame(Yes = 280, No = 120)
  rownames(testdat) <- "B"
  testout <- get_pseudo_ipd_binary(
    binary_agd = testdat,
    format = "unstacked"
  )
  expectout <- data.frame(
    USUBJID = paste0("pseudo_binary_subj_", 1:400),
    ARM = rep("B", 400),
    RESPONSE = c(rep(TRUE, 280), rep(FALSE, 120))
  )

  expect_equal(testout, expectout)
})

test_that("generate pseudo binary IPD given stacked table as input", {
  # example of stacked
  testout <- get_pseudo_ipd_binary(
    binary_agd = data.frame(
      ARM = rep("B", 2),
      RESPONSE = c("YES", "NO"),
      COUNT = c(280, 120)
    ),
    format = "stacked"
  )

  expectout <- data.frame(
    USUBJID = paste0("pseudo_binary_subj_", 1:400),
    ARM = rep("B", 400),
    RESPONSE = c(rep(TRUE, 280), rep(FALSE, 120))
  )

  expect_equal(testout, expectout)
})
