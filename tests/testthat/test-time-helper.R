test_that("setting time factors works", {
  expect_no_condition(
    set_time_conversion(default = "years", days = 1 / 365.25, weeks = 1 / 52.17857, months = 1 / 12, years = 1)
  )
  expect_equal(
    get_time_conversion(),
    c(days = 1 / 365.25, weeks = 1 / 52.17857, months = 1 / 12, years = 1)
  )

  expect_equal(
    get_time_as(1:10, as = NULL),
    1:10
  )

  expect_equal(
    get_time_as(1:10, as = "days"),
    1:10 / (1 / 365.25)
  )

  expect_no_condition(
    set_time_conversion(default = "days", days = 1, weeks = 7, months = 365.25 / 12, years = 365.25)
  )

  expect_equal(
    get_time_conversion(),
    c(days = 1, weeks = 7, months = 365.25 / 12, years = 365.25)
  )

  expect_equal(
    get_time_as(1:10, as = "days"),
    1:10
  )

  expect_equal(
    get_time_as(1:10, as = "weeks"),
    1:10 / 7
  )

  expect_error(get_time_as(letters), "numeric")

  expect_error(
    set_time_conversion(default = "days", days = 0, weeks = Inf, months = NA, years = NaN),
    "Conversion factors must be finite non-zero numerical values: days = 0, weeks = Inf, months = NA, years = NaN"
  )
})
