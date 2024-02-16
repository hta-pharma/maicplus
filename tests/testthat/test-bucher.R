test_that("bucher works as expected", {
  result <- bucher(
    trt = list(est = log(1.1), se = 0.2),
    com = list(est = log(1.3), se = 0.18),
    conf_lv = 0.9
  )

  expected <- list(
    est = -0.1670541,
    se = 0.2690725,
    ci_l = -0.6096389,
    ci_u = 0.2755308,
    pval = 0.5346973
  )
  expect_equal(result, expected)
})
