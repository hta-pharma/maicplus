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
  expect_equal(result, expected)
})
