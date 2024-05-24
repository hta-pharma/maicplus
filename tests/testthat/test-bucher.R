test_that("bucher works as expected", {
  result <- bucher(
    trt = list(est = log(1.1), se = 0.2),
    com = list(est = log(1.3), se = 0.18),
    conf_lv = 0.9
  )

  expected <- list(
    est = -0.1670540846631661247,
    se = 0.26907248094147423467,
    ci_l = -0.60963893085258069604,
    ci_u = 0.27553076152624844664,
    pval = 0.53469725822185854014
  )
  class(expected) <- c("maicplus_bucher", "list")
  expect_equal(result, expected)
})

test_that("find_SE_from_CI works as expected", {
  result <- find_SE_from_CI(CI_lower = 0.55, CI_upper = 0.90, CI_perc = 0.95)
  expected <- 0.12563406495792420192

  expect_equal(result, expected)
})


test_that("bucher print works as expected", {
  bucher_result <- bucher(
    trt = list(est = log(1.1), se = 0.2),
    com = list(est = log(1.3), se = 0.18),
    conf_lv = 0.9
  )
  result <- print(bucher_result, ci_digits = 2, pval_digits = 3)

  expected <- c("-0.17 [-0.61; 0.28]", "0.535")
  names(expected) <- c("result", "pvalue")

  expect_equal(result, expected)
})
