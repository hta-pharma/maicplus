test_that("bucher works as expected", {
  result <- bucher(
    trt = list(est = 1, se = 0.04),
    com = list(est = 1.2, se = 0.08),
    conf_lv = 0.9
  )

  expected <- list(
    est = -0.2,
    se = 0.0894427190999916,
    ci_l = -0.347120180916023,
    ci_u = -0.0528798190839771,
    pval = 0.0253473186774683
  )
  expect_equal(result, expected)
})
