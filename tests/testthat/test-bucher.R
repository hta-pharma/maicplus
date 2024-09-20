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

test_that("bucher works when difference in trt effects is positive", {
  result <- bucher(
    trt = list(est = 1.4, se = 0.09),
    com = list(est = 1.2, se = 0.07),
    conf_lv = 0.9
  )

  est <- 0.2
  se <- sqrt(0.09^2 + 0.07^2)
  q <- qnorm((1 - 0.9) / 2, lower.tail = FALSE)
  pval <- 2 * (1 - stats::pnorm(est, 0, se))
  expected <- list(
    est = est,
    se = se,
    ci_l = est - q * se,
    ci_u = est + q * se,
    pval = pval
  )
  class(expected) <- c("maicplus_bucher", "list")
  expect_equal(result, expected)
})

test_that("bucher errors on bad input", {
  expect_error(bucher(trt = list(est = NA, se = 0.2), com = list(est = 1.3, se = 0.18), conf_lv = 0.9), "est")
  expect_error(bucher(trt = list(est = 1.4, se = NA), com = list(est = 1.3, se = 0.18), conf_lv = 0.9), "se")
  expect_error(bucher(trt = list(est = 1.4, se = 0.2), com = list(est = NA, se = 0.18), conf_lv = 0.9), "est")
  expect_error(bucher(trt = list(est = 1.4, se = 0.2), com = list(est = 1.3, se = NULL), conf_lv = 0.9), "se")
  expect_error(bucher(trt = list(est = 1.4, se = 0.2), com = list(est = 1.3, se = 0.18), conf_lv = 1.9), "conf_lv")
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

  expected <- c("-0.17[-0.61; 0.28]", "0.535")
  names(expected) <- c("result", "pvalue")

  expect_equal(result, expected)
})
