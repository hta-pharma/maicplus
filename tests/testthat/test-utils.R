test_that("set_random_seed works", {
  original_seed <- globalenv()$.Random.seed
  set.seed(123)
  seed_123 <- globalenv()$.Random.seed
  sample(10)
  # Back to a user specified state
  set_random_seed(seed_123)
  expect_equal(seed_123, globalenv()$.Random.seed)
  # Back to the original state
  set_random_seed(original_seed)
  expect_equal(original_seed, globalenv()$.Random.seed)
})

test_that("transform_ratio works as expected", {
  object <- list(est = 2, se = 1, ci_l = 2 - 1.96, ci_u = 2 + 1.96)

  result <- transform_ratio(object)
  expected <- list(
    est = exp(2),
    se = sqrt((exp(1^2) - 1) * exp(2 * 2 + 1^2)),
    ci_l = exp(2 - 1.96),
    ci_u = exp(2 + 1.96)
  )
  expect_equal(result, expected)
})


test_that("transform_absolute works as expected", {
  object <- list(est = 2, se = 1, ci_l = 2 - 1.96, ci_u = 2 + 1.96)

  result <- transform_absolute(object)
  expected <- list(
    est = 200,
    se = 100,
    ci_l = 4,
    ci_u = 396
  )
  expect_equal(result, expected)
})
