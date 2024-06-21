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
